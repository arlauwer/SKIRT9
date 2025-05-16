/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "XRayIonicGasMix.hpp"
#include "ComptonPhaseFunction.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "DipolePhaseFunction.hpp"
#include "FatalError.hpp"
#include "MaterialState.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "Range.hpp"
#include "SnapshotParameter.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"
#include <map>
#include <regex>

////////////////////////////////////////////////////////////////////

// ---- resource data structs ---- //
// these structs are used to store the data read from the resource files
// the final stored data does not necessarily match the data in the resource files

namespace
{
    // All the parameters we want to store, discarded after the setup.
    // These are used to prepare the XRayIonicGasMix.

    struct PhotoAbsorbResource
    {
        PhotoAbsorbResource(const Array& params)
            : Z(params[0]), N(params[1]), n(params[2]), l(params[3]), Eth(params[4]), E0(params[5]), sigma0(params[6]),
              ya(params[7]), P(params[8]), yw(params[9])
        {}
        // resource
        short Z;        // atomic number
        short N;        // number of electrons
        short n;        // principal quantum number of the shell
        short l;        // orbital quantum number of the subshell
        double Eth;     // subshell ionization threshold energy (eV)
        double E0;      // fit parameter (eV)
        double sigma0;  // fit parameter (Mb = 10^-22 m^2)
        double ya;      // fit parameter (1)
        double P;       // fit parameter (1)
        double yw;      // fit parameter (1)

        // calculated
        double Es;
        double sigmamax;

        // pointer
        int ionIndex;  // index of the ion

        // return photo-absorption cross section in m2 for given energy in eV and cross section parameters,
        // without taking into account thermal dispersion
        double photoAbsorbSection(double E) const
        {
            if (E < Eth || E >= 5e5) return 0.;

            double x = E / E0;
            double y = std::sqrt(x * x);
            double xm1 = x - 1.;
            double Q = 5.5 + l - 0.5 * P;
            double F = (xm1 * xm1 + yw * yw) * std::pow(y, -Q) * std::pow(1. + std::sqrt(y / ya), -P);
            return 1e-22 * sigma0 * F;  // from Mb to m2
        }

        // return photo-absorption cross section in m2 for given energy in eV and cross section parameters,
        // approximating thermal dispersion by replacing the steep threshold transition by a sigmoid
        // error function with given parameters (dispersion and maximum value)
        double photoAbsorbThermalSection(double E) const
        {
            if (E <= Eth - 2. * Es) return 0.;
            if (E >= Eth + 2. * Es) return photoAbsorbSection(E);
            return sigmamax * (0.5 + 0.5 * std::erf((E - Eth) / Es));
        }

        // calculate the parameters for the sigmoid function approximating the convolution with a Gaussian
        // at the threshold energy for each cross section record, and store the result into a temporary vector;
        // the information includes the thermal energy dispersion at the threshold energy and
        // the intrinsic cross section at the threshold energy plus twice this energy dispersion
        void setThermalDispersion(double vtherm)
        {
            Es = Eth * vtherm / Constants::c();
            sigmamax = photoAbsorbSection(Eth + 2. * Es);
        }
    };

    struct FluorescenceResource
    {
        FluorescenceResource(const Array& params)
            : Z(params[0]), N(params[1]), n(params[2]), l(params[3]), omega(params[4]), E(params[5]), W(params[6])
        {}

        // resource
        short Z;       // atomic number
        short N;       // number of electrons
        short n;       // principal quantum number of the corresponding photo-absorption
        short l;       // orbital quantum number of the corresponding photo-absorption
        double omega;  // fluorescence yield (1)
        double E;      // (central) energy of the emitted photon (eV)
        double W;      // FWHM of the Lorentz shape for the emitted photon (eV), or zero

        // pointer
        int papIndex;  // index of the corresponding photo-absorption
    };

    struct BoundTransitionResource
    {
        BoundTransitionResource(const Array& params)
            : Z(params[0]), N(params[1]), res_upperIndex(params[2]), res_lowerIndex(params[3]), lam(params[4]),
              A(params[5]), gul(params[6])
        {}

        // resource
        short Z;             // atomic number
        short N;             // number of electrons
        int res_upperIndex;  // index of the upper energy level in the Chianti order (resource)
        int res_lowerIndex;  // index of the lower energy level in the Chianti order (resource)
        double lam;          // wavelength of the transition (m)
        double A;            // Einstein A coefficient (s^-1)
        double gul;          // gu/gl

        // calculated
        int upperIndex;  // index of the upper energy level in the BoundLevel order
        int lowerIndex;  // index of the lower energy level in the BoundLevel order

        // pointer
        int ionIndex;  // index of the ion
    };

    struct BoundLevelInput
    {
        BoundLevelInput(const Array& params) : Z(params[0]), N(params[1]), index(params[2]), pop(params[3]) {}

        // user input
        short Z;     // atomic number
        short N;     // number of electrons
        int index;   // Chianti index of the level
        double pop;  // level population
    };
}

////////////////////////////////////////////////////////////////////

namespace
{
    // ---- common helper functions ---- //

    static constexpr double defaultTemperature = 1;

    // map of element names to atomic numbers
    const std::map<std::string, short> elementToZ = {
        {"H", 1},   {"He", 2},  {"Li", 3},  {"Be", 4},  {"B", 5},   {"C", 6},   {"N", 7},  {"O", 8},
        {"F", 9},   {"Ne", 10}, {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16},
        {"Cl", 17}, {"Ar", 18}, {"K", 19},  {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24},
        {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30}};

    // masses of the elements in amu up to Z=30
    const vector<double> elementMasses = {1.0079,  4.0026, 6.941,   9.01218, 10.81,   12.011,  14.0067, 15.9994,
                                          18.9984, 20.179, 22.9898, 24.305,  26.9815, 28.0855, 30.9738, 32.06,
                                          35.453,  39.948, 39.0983, 40.08,   44.9559, 47.9,    50.9415, 51.996,
                                          54.938,  55.847, 58.9332, 58.7,    63.546,  65.38};

    // map ionString (eg. "Fe+4") to atomic number Z and number of electrons N
    XRayIonicGasMix::Ion stringToIon(string ionString)
    {
        ionString = StringUtils::squeeze(ionString);
        // split ion string into element symbol and ionization number
        std::regex pattern("([A-Za-z]+)\\+?([0-9]*)");
        std::smatch match;
        if (!std::regex_match(ionString, match, pattern)) throw FATALERROR("Invalid ion format: " + ionString);

        // get ion parameters
        int Z = elementToZ.at(match[1].str());
        int N = Z - std::stoi(match[2].str());
        if (N < 1 || N > Z) throw FATALERROR("Invalid ionization stage: " + ionString);
        if (Z < 1 || Z > 30) throw FATALERROR("Invalid atomic number: " + ionString);

        XRayIonicGasMix::Ion ion;
        ion.Z = Z;
        ion.N = N;
        return ion;
    }

    // convert photon energy in eV to and from wavelength in m (same conversion in both directions)
    constexpr double wavelengthToFromEnergy(double x)
    {
        constexpr double front = Constants::h() * Constants::c() / Constants::Qelectron();
        return front / x;
    }

    // multiplicator to convert energy in keV to scaled energy E / (m_e c^2)
    constexpr double keVtoScaledEnergy =
        (1e3 * Constants::Qelectron()) / (Constants::Melectron() * Constants::c() * Constants::c());

    // multiplicator to convert scaled energy to energy in units of 12.4 keV
    constexpr double scaledEnergyTo12keV =
        (Constants::Melectron() * Constants::c() * Constants::c()) / (12.4e3 * Constants::Qelectron());

    // convert wavelength to scaled photon energy: h nu / (m_e c^2)
    constexpr double scaledEnergy(double lambda)
    {
        constexpr double front = Constants::h() / Constants::Melectron() / Constants::c();
        return front / lambda;
    }

    // ---- hardcoded configuration constants ---- //

    // wavelength range over which our cross sections may be nonzero
    constexpr Range nonZeroRange(wavelengthToFromEnergy(500e3), wavelengthToFromEnergy(4.3));

    // number of wavelengths per dex in high-resolution grid
    constexpr size_t numWavelengthsPerDex = 2500;

    // discretization of the phase function over scattering angle: theta from 0 to pi, index t
    constexpr size_t numTheta = 361;
    constexpr size_t maxTheta = numTheta - 1;
    constexpr double deltaTheta = M_PI / maxTheta;

    // load data from resource file with N columns into a vector of structs of type S that can be constructed
    // from an array with N elements, and return that vector
    template<class S, int N>
    vector<S> loadStruct(const SimulationItem* item, string filename, const vector<XRayIonicGasMix::Ion>& ions,
                         bool resource = true)
    {
        if (filename.empty()) return vector<S>();

        vector<S> result;
        TextInFile infile(item, filename, "", resource);
        for (int i = 0; i != N; ++i) infile.addColumn(string());
        Array row;
        while (infile.readRow(row))
        {
            S s(row);
            for (const auto& ion : ions)
            {
                if (ion.Z == s.Z && ion.N == s.N)
                {
                    result.push_back(s);
                    break;
                }
            }
        }
        return result;
    }

    // ---- bound-electron scattering resources ---- //

    // load data from resource file with N columns into a vector of N arrays, and return that vector;
    // each of the arrays is resized to remove trailing NaN values, if applicable
    vector<Array> loadColumns(int N, const SimulationItem* item, string filename, string description)
    {
        TextInFile infile(item, filename, description, false);
        for (int i = 0; i != N; ++i) infile.addColumn(string());
        vector<Array> columns = infile.readAllColumns();

        // clip any columns with trailing NaNs
        for (Array& column : columns)
        {
            size_t n = column.size();
            while (n && std::isnan(column[n - 1])) --n;
            if (n != column.size())
            {
                // we need to make a copy because resizing an array clears its contents
                Array copy = column;
                column.resize(n);
                for (size_t i = 0; i != n; ++i) column[i] = copy[i];
            }
        }
        return columns;
    }

    double gaussian(double x, double mu, double stddev)
    {
        constexpr double front = 0.25 * M_SQRT2 * M_2_SQRTPI;
        double u = (x - mu) / stddev;
        return front * exp(-0.5 * u * u);
    }
}

////////////////////////////////////////////////////////////////////

// ---- base class for scattering helpers ---- //

class XRayIonicGasMix::ScatteringHelper
{
public:
    virtual ~ScatteringHelper() {}

    // return scattering cross section for atom in m2
    virtual double sectionSca(double lambda, int Z) const = 0;

    // peel-off unpolarized scattering event: override this in helpers that don't support polarization
    virtual void peeloffScattering(double& /*I*/, double& /*lambda*/, int /*Z*/, Direction /*bfk*/,
                                   Direction /*bfkobs*/) const
    {
        // default implementation does nothing
    }

    // perform unpolarized scattering event: override this in helpers that don't support polarization
    virtual Direction performScattering(double& /*lambda*/, int /*Z*/, Direction /*bfk*/) const
    {
        // default implementation returns null vector
        return Direction();
    }

    // peel-off polarized scattering event: override this in helpers that do support polarization
    virtual void peeloffScattering(double& I, double& /*Q*/, double& /*U*/, double& /*V*/, double& lambda, int Z,
                                   Direction bfk, Direction bfkobs, Direction /*bfky*/,
                                   const StokesVector* /*sv*/) const
    {
        // default implementation calls unpolarized version
        peeloffScattering(I, lambda, Z, bfk, bfkobs);
    }

    // perform polarized scattering event: override this in helpers that do support polarization
    virtual Direction performScattering(double& lambda, int Z, Direction bfk, StokesVector* /*sv*/) const
    {
        // default implementation calls unpolarized version
        return performScattering(lambda, Z, bfk);
    }
};

////////////////////////////////////////////////////////////////////

// ---- no scattering helper ---- //

namespace
{
    // this helper does nothing; it is used as a stub in case there is no scattering of a given type
    class NoScatteringHelper : public XRayIonicGasMix::ScatteringHelper
    {
    public:
        NoScatteringHelper(SimulationItem* /*item*/) {}

        double sectionSca(double /*lambda*/, int /*Z*/) const override { return 0.; }
    };
}

////////////////////////////////////////////////////////////////////

// ---- free-electron Compton scattering helper ---- //

namespace
{
    // transition wavelength from Compton to Thomson scattering
    constexpr double comptonWL = wavelengthToFromEnergy(100.);  // 0.1 keV or 12.4 nm

    // this helper forwards all calls to an external helper class for regular Compton scattering
    // (or Thomson scattering for lower energies, because Compton becomes numerically unstable)
    class FreeComptonHelper : public XRayIonicGasMix::ScatteringHelper
    {
    private:
        ComptonPhaseFunction _cpf;
        DipolePhaseFunction _dpf;

    public:
        FreeComptonHelper(SimulationItem* item)
        {
            auto random = item->find<Random>();
            _cpf.initialize(random);
            _dpf.initialize(random);
        }

        double sectionSca(double lambda, int Z) const override
        {
            double sigma = Z * Constants::sigmaThomson();
            if (lambda < comptonWL) sigma *= _cpf.sectionSca(lambda);
            return sigma;
        }

        void peeloffScattering(double& I, double& lambda, int /*Z*/, Direction bfk, Direction bfkobs) const override
        {
            if (lambda < comptonWL)
            {
                double Q, U, V;
                _cpf.peeloffScattering(I, Q, U, V, lambda, bfk, bfkobs, Direction(), nullptr);
            }
            else
            {
                double Q, U, V;
                _dpf.peeloffScattering(I, Q, U, V, bfk, bfkobs, Direction(), nullptr);
            }
        }

        Direction performScattering(double& lambda, int /*Z*/, Direction bfk) const override
        {
            return lambda < comptonWL ? _cpf.performScattering(lambda, bfk, nullptr)
                                      : _dpf.performScattering(bfk, nullptr);
        }
    };
}

////////////////////////////////////////////////////////////////////

// ---- free-electron Compton with polarization scattering helper ---- //

namespace
{
    // this helper forwards all calls to an external helper class for Compton scattering
    // (or Thomson scattering for lower energies) with support for polarization
    class FreeComptonWithPolarizationHelper : public XRayIonicGasMix::ScatteringHelper
    {
    private:
        ComptonPhaseFunction _cpf;
        DipolePhaseFunction _dpf;

    public:
        FreeComptonWithPolarizationHelper(SimulationItem* item)
        {
            auto random = item->find<Random>();
            _cpf.initialize(random, true);
            _dpf.initialize(random, true);
        }

        double sectionSca(double lambda, int Z) const override
        {
            double sigma = Z * Constants::sigmaThomson();
            if (lambda < comptonWL) sigma *= _cpf.sectionSca(lambda);
            return sigma;
        }

        void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, int /*Z*/, Direction bfk,
                               Direction bfkobs, Direction bfky, const StokesVector* sv) const override
        {
            lambda < comptonWL ? _cpf.peeloffScattering(I, Q, U, V, lambda, bfk, bfkobs, bfky, sv)
                               : _dpf.peeloffScattering(I, Q, U, V, bfk, bfkobs, bfky, sv);
        }

        Direction performScattering(double& lambda, int /*Z*/, Direction bfk, StokesVector* sv) const override
        {
            return lambda < comptonWL ? _cpf.performScattering(lambda, bfk, sv) : _dpf.performScattering(bfk, sv);
        }
    };
}

////////////////////////////////////////////////////////////////////

// ---- bound-electron Compton scattering helper ---- //

namespace
{
    // returns the inverse Compton factor for a given scaled energy and scattering angle cosine
    constexpr double inverseComptonFactor(double x, double costheta)
    {
        return 1. + x * (1. - costheta);
    }

    // returns the Compton factor for a given scaled energy and scattering angle cosine
    constexpr double comptonFactor(double x, double costheta)
    {
        return 1. / inverseComptonFactor(x, costheta);
    }

    // returns the value interpolated from the specified table as a function of the momentum transfer parameter
    // q = (E/12.4 keV) sin(theta/2), given the scaled energy x and the sine;
    // logarithmic interpolation is used except for q values near zero
    double interpolateQ(double x, double sintheta2, const Array& qv, const Array& fv)
    {
        double q = scaledEnergyTo12keV * x * sintheta2;
        if (q < 1e-3) return NR::clampedValue<NR::interpolateLinLin>(q, qv, fv);
        return NR::clampedValue<NR::interpolateLogLog>(q, qv, fv);
    }

    // this helper implements bound-electron Compton scattering
    class BoundComptonHelper : public XRayIonicGasMix::ScatteringHelper
    {
    private:
        // resources loaded from file
        vector<Array> _CSv;  // 0: E (keV->1); 1-30: bound Compton cross sections (cm2->m2)
        vector<Array> _SFv;  // 0: q (1); 1-30: incoherent scattering functions (1)
        vector<Array> _CPv;  // 0: E (keV->1); 1-30: pdf for target electron momentum (1)
        vector<Array> _IBv;  // 0: E (keV->1) ionisation energy of the outer subshell electrons

        // precalculated cumulative distributions for target electron momentum
        Range _cumRange;
        vector<Array> _cumCPv;  // 0: E axis; 1-30: cumulative pdf for target electron momentum

        // precalculated discretizations
        Array _costhetav = Array(numTheta);
        Array _sinthetav = Array(numTheta);
        Array _sin2thetav = Array(numTheta);
        Array _sintheta2v = Array(numTheta);

        // cache
        Random* _random{nullptr};

    public:
        BoundComptonHelper(SimulationItem* item)
        {
            // load bound Compton cross sections
            _CSv = loadColumns(30 + 1, item, "XRay_CS.txt", "bound Compton data");
            _CSv[0] *= keVtoScaledEnergy;                      // convert from keV to 1
            for (size_t Z = 1; Z <= 30; ++Z) _CSv[Z] *= 1e-4;  // convert from cm2 to m2

            // load incoherent scattering functions
            _SFv = loadColumns(30 + 1, item, "XRay_SF.txt", "bound Compton data");

            // load pdfs for projected momentum of target electron
            _CPv = loadColumns(30 + 1, item, "XRay_CP.txt", "bound Compton data");
            _CPv[0] *= keVtoScaledEnergy;  // convert from keV to 1

            // load ionization energies
            _IBv = loadColumns(1, item, "XRay_IB.txt", "bound Compton data");
            _IBv[0] *= keVtoScaledEnergy;  // convert from keV to 1

            // precalculate cumulative distributions for target electron momentum
            _cumRange.set(_CPv[0][0], _CPv[0][_CPv[0].size() - 1]);
            Array xv, pv, Pv;
            NR::cdf<NR::interpolateLinLin>(xv, pv, Pv, _CPv[0], _CPv[1], _cumRange);
            _cumCPv.push_back(xv);
            _cumCPv.push_back(Pv);
            for (size_t Z = 2; Z <= 30; ++Z)
            {
                NR::cdf<NR::interpolateLinLin>(xv, pv, Pv, _CPv[0], _CPv[1], _cumRange);
                _cumCPv.push_back(Pv);
            }

            // construct a theta grid and precalculate values used in generateCosineFromPhaseFunction()
            // to accelerate construction of the cumulative phase function distribution
            for (size_t t = 0; t != numTheta; ++t)
            {
                double theta = t * deltaTheta;
                _costhetav[t] = cos(theta);
                _sinthetav[t] = sin(theta);
                _sin2thetav[t] = _sinthetav[t] * _sinthetav[t];
                _sintheta2v[t] = sin(0.5 * theta);
            }

            // cache random nr generator
            _random = item->find<Random>();
        }

        double sectionSca(double lambda, int Z) const override
        {
            // interpolate from table, and:
            // - below lower table limit: cross section must be zero so don't clamp values
            // - above upper table limit: does not matter because this limit coincides with the global upper limit
            return NR::value<NR::interpolateLogLog>(scaledEnergy(lambda), _CSv[0], _CSv[Z]);
        }

    private:
        double phaseFunctionValue(double x, double costheta, int Z) const
        {
            constexpr double norm = 3. / 4. * Constants::sigmaThomson();
            double C = comptonFactor(x, costheta);
            double sin2theta = (1 - costheta) * (1 + costheta);
            double phase = C * C * C + C - C * C * sin2theta;
            double section = NR::value<NR::interpolateLogLog>(x, _CSv[0], _CSv[Z]);
            double sintheta2 = sqrt(0.5 * (1 - costheta));
            double incoherent = interpolateQ(x, sintheta2, _SFv[0], _SFv[Z]);
            return norm / section * phase * incoherent;
        }

        double generateCosineFromPhaseFunction(double x, double Z) const
        {
            // construct the normalized cumulative phase function distribution for this x
            Array thetaXv;
            NR::cdf(thetaXv, maxTheta, [this, x, Z](int t) {
                t += 1;
                double C = comptonFactor(x, _costhetav[t]);
                double phase = C * C * C + C - C * C * _sin2thetav[t];
                double incoherent = interpolateQ(x, _sintheta2v[t], _SFv[0], _SFv[Z]);
                return phase * incoherent * _sinthetav[t];
            });

            // draw a random cosine from this distribution
            return _random->cdfLinLin(_costhetav, thetaXv);
        }

        // sample a target electron momentum from the distribution with the given maximum
        double sampleMomentum(double pmax, double Z) const
        {
            // maximum momentum is below the range of the tabulated pdf -> simply return the maximum momentum
            // (we estimate that this happens for less than 0.1 % of the events)
            if (pmax <= _cumRange.min()) return pmax;

            // maximum momentum is on the left side of the peak in the tabulated pdf;
            // using the rejection technique on the full-range pdf is very inefficient
            // because the majority of the generated samples would be rejected
            // --> reconstruct a cumulative pdf with the appropriate range and use numerical inversion
            // (we estimate that this happens for less than 10% of the events)
            if (pmax <= _cumRange.mid())
            {
                Array xv, pv, Pv;
                NR::cdf<NR::interpolateLinLin>(xv, pv, Pv, _CPv[0], _CPv[Z], Range(_cumRange.min(), pmax));
                return _random->cdfLinLin(xv, Pv);
            }

            // maximum momentum is on the right side of the peak in the tabulated pdf, possibly even out of range;
            // using the rejection technique on top of numerical inversion for the full-range pdf now is efficient
            // and quite fast because we can use the precalculated cumulative pdf
            while (true)
            {
                double p = _random->cdfLinLin(_cumCPv[0], _cumCPv[Z]);
                if (p <= pmax) return p;
            }
        }

        // returns the augmented inverse Compton factor
        double augmentedInverseComptonFactor(double x, double costheta, double Z) const
        {
            // precalculate some values
            double costheta1 = 1. - costheta;
            double sintheta22 = 2. * sqrt(0.5 * costheta1);  // twice the half-angle sine

            // calculate the maximum target electron momentum (in scaled energy units)
            double b = _IBv[0][Z - 1];  // scaled ionization energy
            double xminb = (x - b);
            double pmax = (x * xminb * costheta1 - b) / (xminb * sintheta22);

            // sample a target electron momentum from the distribution with the given maximum
            double p = sampleMomentum(pmax, Z);

            // calculate the augmented inverse Compton factor
            return 1. + x * costheta1 - p * sintheta22;
        }

    public:
        void peeloffScattering(double& I, double& lambda, int Z, Direction bfk, Direction bfkobs) const override
        {
            double x = scaledEnergy(lambda);

            // calculate the value of the phase function
            double costheta = Vec::dot(bfk, bfkobs);
            double value = phaseFunctionValue(x, costheta, Z);

            // accumulate the weighted sum in the intensity
            I += value;

            // adjust the wavelength
            lambda *= augmentedInverseComptonFactor(x, costheta, Z);
        }

        Direction performScattering(double& lambda, int Z, Direction bfk) const override
        {
            double x = scaledEnergy(lambda);

            // sample a scattering angle from the phase function
            double costheta = generateCosineFromPhaseFunction(x, Z);

            // adjust the wavelength
            lambda *= augmentedInverseComptonFactor(x, costheta, Z);

            // determine the new propagation direction
            return _random->direction(bfk, costheta);
        }
    };
}

////////////////////////////////////////////////////////////////////

// ---- smooth Rayleigh scattering helper ---- //

namespace
{
    // this helper implements smooth Rayleigh scattering;
    // below the energy limit of the tabulated data, use Thomson scattering instead
    class SmoothRayleighHelper : public XRayIonicGasMix::ScatteringHelper
    {
    private:
        vector<Array> _RSSv;  // 0: E (keV->1); 1-30: smooth Rayleigh cross sections (cm2->m2)
        vector<Array> _FFv;   // 0: q (1); 1-30: atomic form factors (1)
        Random* _random{nullptr};
        DipolePhaseFunction _dpf;

        // precalculated discretizations
        Array _costhetav = Array(numTheta);
        Array _cos2thetav = Array(numTheta);
        Array _sinthetav = Array(numTheta);
        Array _sintheta2v = Array(numTheta);

    public:
        SmoothRayleighHelper(SimulationItem* item)
        {
            // load smooth Rayleigh cross sections
            _RSSv = loadColumns(30 + 1, item, "XRay_RSS.txt", "smooth Rayleigh data");
            _RSSv[0] *= keVtoScaledEnergy;                   // convert from keV to 1
            for (int Z = 1; Z <= 30; ++Z) _RSSv[Z] *= 1e-4;  // convert from cm2 to m2

            // load atomic form factors
            _FFv = loadColumns(30 + 1, item, "XRay_FF.txt", "smooth Rayleigh data");

            // cache random nr generator and initialize the Thomson helper
            _random = item->find<Random>();
            _dpf.initialize(_random);

            // construct a theta grid and precalculate values used in generateCosineFromPhaseFunction()
            // to accelerate construction of the cumulative phase function distribution
            for (size_t t = 0; t != numTheta; ++t)
            {
                double theta = t * deltaTheta;
                _costhetav[t] = cos(theta);
                _cos2thetav[t] = _costhetav[t] * _costhetav[t];
                _sinthetav[t] = sin(theta);
                _sintheta2v[t] = sin(0.5 * theta);
            }
        }

        double sectionSca(double lambda, int Z) const override
        {
            // interpolate from table, and:
            // - below lower table limit: use Z^2 * Thomson scattering
            // - above upper table limit: does not matter because this limit coincides with the global upper limit
            double x = scaledEnergy(lambda);
            if (x < _RSSv[0][0]) return Z * Z * Constants::sigmaThomson();
            return NR::value<NR::interpolateLogLog>(x, _RSSv[0], _RSSv[Z]);
        }

    private:
        double phaseFunctionValue(double x, double costheta, int Z) const
        {
            constexpr double norm = 3. / 4. * Constants::sigmaThomson();
            double phase = 1. + costheta * costheta;
            double section = NR::value<NR::interpolateLogLog>(x, _RSSv[0], _RSSv[Z]);
            double sintheta2 = sqrt(0.5 * (1 - costheta));
            double form = interpolateQ(x, sintheta2, _FFv[0], _FFv[Z]);
            return norm / section * phase * form * form;
        }

        double generateCosineFromPhaseFunction(double x, double Z) const
        {
            // construct the normalized cumulative phase function distribution for this x
            Array thetaXv;
            NR::cdf(thetaXv, maxTheta, [this, x, Z](int t) {
                t += 1;
                double phase = 1. + _cos2thetav[t];
                double form = interpolateQ(x, _sintheta2v[t], _FFv[0], _FFv[Z]);
                return phase * form * form * _sinthetav[t];
            });

            // draw a random cosine from this distribution
            return _random->cdfLinLin(_costhetav, thetaXv);
        }

    public:
        void peeloffScattering(double& I, double& lambda, int Z, Direction bfk, Direction bfkobs) const override
        {
            double x = scaledEnergy(lambda);

            // for low energies use Thomson scattering
            if (x < _RSSv[0][0])
            {
                double Q, U, V;
                _dpf.peeloffScattering(I, Q, U, V, bfk, bfkobs, Direction(), nullptr);
            }

            // otherwise use Rayleigh scattering
            else
            {
                // calculate the value of the phase function
                double costheta = Vec::dot(bfk, bfkobs);
                double value = phaseFunctionValue(x, costheta, Z);

                // accumulate the weighted sum in the intensity
                I += value;
            }
        }

        Direction performScattering(double& lambda, int Z, Direction bfk) const override
        {
            double x = scaledEnergy(lambda);

            // for low energies use Thomson scattering
            if (x < _RSSv[0][0]) return _dpf.performScattering(bfk, nullptr);

            // otherwise use Rayleigh scattering
            return _random->direction(bfk, generateCosineFromPhaseFunction(x, Z));
        }
    };
}

////////////////////////////////////////////////////////////////////

// ---- anomalous Rayleigh scattering helper ---- //

namespace
{
    // this helper implements anomalous Rayleigh scattering
    // below the energy limit of the tabulated data, use Thomson scattering instead
    class AnomalousRayleighHelper : public XRayIonicGasMix::ScatteringHelper
    {
    private:
        vector<Array> _RSAv;  // 2*Z: E (keV->1); 2*Z+1: anomalous Rayleigh cross sections (cm2->m2)
        vector<Array> _FFv;   // 0: q (1); 1-30: atomic form factors (1)
        vector<Array> _F1v;   // 2*Z: E (keV->1); 2*Z+1: Real anomalous scattering function (1)
        vector<Array> _F2v;   // 2*Z: E (keV->1); 2*Z+1: Imaginary anomalous scattering function (1)
        Random* _random{nullptr};
        DipolePhaseFunction _dpf;

        // precalculated discretizations
        Array _costhetav = Array(numTheta);
        Array _cos2thetav = Array(numTheta);
        Array _sinthetav = Array(numTheta);
        Array _sintheta2v = Array(numTheta);

    public:
        AnomalousRayleighHelper(SimulationItem* item)
        {
            // load anomalous Rayleigh cross sections, atomic form factors and anomalous scattering functions
            _RSAv = loadColumns(2 * 30 + 2, item, "XRay_RSA.txt", "anomalous Rayleigh data");
            _FFv = loadColumns(30 + 1, item, "XRay_FF.txt", "anomalous Rayleigh data");
            _F1v = loadColumns(2 * 30 + 2, item, "XRay_F1.txt", "anomalous Rayleigh data");
            _F2v = loadColumns(2 * 30 + 2, item, "XRay_F2.txt", "anomalous Rayleigh data");

            // convert units
            for (size_t Z = 1; Z <= 30; ++Z)
            {
                _RSAv[2 * Z] *= keVtoScaledEnergy;  // convert from keV to 1
                _F1v[2 * Z] *= keVtoScaledEnergy;   // convert from keV to 1
                _F2v[2 * Z] *= keVtoScaledEnergy;   // convert from keV to 1
                _RSAv[2 * Z + 1] *= 1e-4;           // convert from cm2 to m2
            }

            // cache random nr generator and initialize the Thomson helper
            _random = item->find<Random>();
            _dpf.initialize(_random);

            // construct a theta grid and precalculate values used in generateCosineFromPhaseFunction()
            // to accelerate construction of the cumulative phase function distribution
            for (size_t t = 0; t != numTheta; ++t)
            {
                double theta = t * deltaTheta;
                _costhetav[t] = cos(theta);
                _cos2thetav[t] = _costhetav[t] * _costhetav[t];
                _sinthetav[t] = sin(theta);
                _sintheta2v[t] = sin(0.5 * theta);
            }
        }

        double sectionSca(double lambda, int Z) const override
        {
            // interpolate from table, and:
            // - below lower table limit: use Z^2 * Thomson scattering
            // - above upper table limit: use clamped value
            double x = scaledEnergy(lambda);
            if (x < _RSAv[2 * Z][0]) return Z * Z * Constants::sigmaThomson();
            return NR::clampedValue<NR::interpolateLogLog>(x, _RSAv[2 * Z], _RSAv[2 * Z + 1]);
        }

    private:
        double phaseFunctionValue(double x, double costheta, int Z) const
        {
            constexpr double norm = 3. / 4. * Constants::sigmaThomson();
            double phase = 1. + costheta * costheta;
            double section = NR::clampedValue<NR::interpolateLogLog>(x, _RSAv[2 * Z], _RSAv[2 * Z + 1]);
            double sintheta2 = sqrt(0.5 * (1 - costheta));
            double form = interpolateQ(x, sintheta2, _FFv[0], _FFv[Z]);
            double form1 = NR::clampedValue<NR::interpolateLogLin>(x, _F1v[2 * Z], _F1v[2 * Z + 1]);  // negative values
            double form2 = NR::clampedValue<NR::interpolateLogLog>(x, _F2v[2 * Z], _F2v[2 * Z + 1]);
            double formsum = form + form1;
            return norm / section * phase * (formsum * formsum + form2 * form2);
        }

        double generateCosineFromPhaseFunction(double x, double Z) const
        {
            // construct the normalized cumulative phase function distribution for this x
            Array thetaXv;
            NR::cdf(thetaXv, maxTheta, [this, x, Z](int t) {
                t += 1;
                double phase = 1. + _cos2thetav[t];
                double form = interpolateQ(x, _sintheta2v[t], _FFv[0], _FFv[Z]);
                double form1 = NR::clampedValue<NR::interpolateLogLin>(x, _F1v[2 * Z], _F1v[2 * Z + 1]);
                double form2 = NR::clampedValue<NR::interpolateLogLog>(x, _F2v[2 * Z], _F2v[2 * Z + 1]);
                double formsum = form + form1;
                return phase * (formsum * formsum + form2 * form2) * _sinthetav[t];
            });

            // draw a random cosine from this distribution
            return _random->cdfLinLin(_costhetav, thetaXv);
        }

    public:
        void peeloffScattering(double& I, double& lambda, int Z, Direction bfk, Direction bfkobs) const override
        {
            double x = scaledEnergy(lambda);

            // for low energies use Thomson scattering
            if (x < _RSAv[2 * Z][0])
            {
                double Q, U, V;
                _dpf.peeloffScattering(I, Q, U, V, bfk, bfkobs, Direction(), nullptr);
            }

            // otherwise use Rayleigh scattering
            else
            {
                // calculate the value of the phase function
                double costheta = Vec::dot(bfk, bfkobs);
                double value = phaseFunctionValue(x, costheta, Z);

                // accumulate the weighted sum in the intensity
                I += value;
            }
        }

        Direction performScattering(double& lambda, int Z, Direction bfk) const override
        {
            double x = scaledEnergy(lambda);

            // for low energies use Thomson scattering
            if (x < _RSAv[2 * Z][0]) return _dpf.performScattering(bfk, nullptr);

            // otherwise use Rayleigh scattering
            return _random->direction(bfk, generateCosineFromPhaseFunction(x, Z));
        }
    };
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMix::setupSelfBefore()
{
    MaterialMix::setupSelfBefore();
    auto config = find<Configuration>();

    // ---- read ions ---- //
    string ionNamesStr = StringUtils::squeeze(ionNames());
    if (ionNamesStr.empty()) throw FATALERROR("No ions specified for XRayIonicGasMix");
    auto ionNameList = StringUtils::split(ionNamesStr, ",");

    for (size_t i = 0; i < ionNameList.size(); ++i)
    {
        string ionName = ionNameList[i];
        auto ion = stringToIon(ionName);

        ion.mass = elementMasses[ion.Z - 1];
        ion.vth = vtherm(ion.mass);
        ion.abund = abundances()[i];  // store for convenience

        _ions.push_back(ion);
    }
    _numIons = _ions.size();

    // ---- load resources ---- //

    // resources
    auto pa_resources = loadStruct<PhotoAbsorbResource, 10>(this, "PA_data.txt", _ions, false);
    auto fl_resources = loadStruct<FluorescenceResource, 7>(this, "FL_data.txt", _ions, false);
    vector<BoundTransitionResource> bt_resources;
    vector<BoundLevelInput> bl_input;
    if (levelPopulations())
    {
        bt_resources = loadStruct<BoundTransitionResource, 7>(this, "BB_data.txt", _ions, false);
        // user input
        bl_input = loadStruct<BoundLevelInput, 4>(this, levelPopulationsFilename(), _ions, false);
    }

    // ---- read bound-bound levels ---- //
    for (auto& blr : bl_input)
    {
        BoundLevel bl;
        bl.pop = blr.pop;
        _boundLevels.push_back(bl);
    }
    _numBLevels = _boundLevels.size();

    // ---- for each ion ---- //
    for (int i = 0; i < _numIons; i++)
    {
        auto& ion = _ions[i];

        // ---- for each photo-absorption ---- //
        int numPa = pa_resources.size();
        for (int p = 0; p < numPa; p++)
        {
            PhotoAbsorbResource& pa_r = pa_resources[p];

            if (pa_r.Z == ion.Z && pa_r.N == ion.N)
            {
                pa_r.setThermalDispersion(vtherm(ion.mass));

                pa_r.ionIndex = i;

                // ---- for each fluorescence ---- //
                for (FluorescenceResource& fl_r : fl_resources)
                {
                    if (fl_r.Z == pa_r.Z && fl_r.N == pa_r.N && fl_r.n == pa_r.n && fl_r.l == pa_r.l)
                    {
                        Fluorescence fl;
                        fl.E = fl_r.E;
                        fl.W = fl_r.W;
                        fl.ionIndex = i;
                        fl_r.papIndex = p;

                        // if no width then convert to wavelengthF
                        if (!fl.W)
                            fl.E = wavelengthToFromEnergy(fl.E);
                        else
                            fl.W /= 2.;

                        _fluorescence.push_back(fl);
                    }
                }
            }
        }

        // ---- for each bound-bound transition ---- //
        for (auto& bt_r : bt_resources)
        {
            if (bt_r.Z == ion.Z && bt_r.N == ion.N)
            {
                bt_r.lam = 1.23984193e-6 / bt_r.lam;  // convert from eV to m

                bt_r.ionIndex = i;

                BoundTransition bt;
                bt.lam = bt_r.lam;
                bt.A = bt_r.A;
                bt.gul = bt_r.gul;
                bt.ionIndex = i;
                for (int b = 0; b < _numBLevels; b++)
                {
                    const BoundLevelInput& bl_in = bl_input[b];
                    // convert the upper and lower indices from Chianti index to the vector index in _boundLevels
                    if (bl_in.Z == bt_r.Z && bl_in.N == bt_r.N && bl_in.index == bt_r.res_upperIndex)
                    {
                        bt_r.upperIndex = b;
                        bt.upperIndex = b;
                    }

                    if (bl_in.Z == bt_r.Z && bl_in.N == bt_r.N && bl_in.index == bt_r.res_lowerIndex)
                    {
                        bt_r.lowerIndex = b;
                        bt.lowerIndex = b;
                    }
                }

                _boundTransitions.push_back(bt);
            }
        }
    }

    _numFluo = _fluorescence.size();
    _numBTransitions = _boundTransitions.size();

    // create scattering helpers depending on the user-configured implementation type;
    // the respective helper constructors load the required bound-electron scattering resources
    switch (scatterBoundElectrons())
    {
        case BoundElectrons::None:
            _ray = new NoScatteringHelper(this);
            _com = new NoScatteringHelper(this);
            break;
        case BoundElectrons::Free:
            _ray = new NoScatteringHelper(this);
            _com = new FreeComptonHelper(this);
            break;
        case BoundElectrons::FreeWithPolarization:
            _ray = new NoScatteringHelper(this);
            _com = new FreeComptonWithPolarizationHelper(this);
            break;
        case BoundElectrons::Good:
            _ray = new SmoothRayleighHelper(this);
            _com = new BoundComptonHelper(this);
            break;
        case BoundElectrons::Exact:
            _ray = new AnomalousRayleighHelper(this);
            _com = new BoundComptonHelper(this);
            break;
    }

    // ---- wavelength grid ---- //

    // construct a wavelength grid for sampling cross sections containing a merged set of grid points
    // in the relevant wavelength range (intersection of simulation range and nonzero range):
    //  - a fine grid in log space that provides sufficient resolution for most applications
    //  - all specific wavelengths mentioned in the configuration of the simulation (grids, normalizations, ...)
    //    ensuring that the cross sections are calculated at exactly these wavelengths
    //  - 7 extra wavelength points around the threshold energies for all transitions,
    //    placed at -2, -4/3, -2/3, 0, 2/3, 4/3, 2 times the thermal energy dispersion

    // we first gather all the wavelength points, in arbitrary order, and then sort them
    vector<double> lambdav;
    lambdav.reserve(5 * numWavelengthsPerDex);

    // get the relevant range (intersection of simulation range and nonzero range)
    Range range = config->simulationWavelengthRange();
    range.intersect(nonZeroRange);

    // add a fine grid in log space;
    // use integer multiples as logarithmic grid points so that the grid is stable for changing wavelength ranges
    constexpr double numPerDex = numWavelengthsPerDex;  // converted to double to avoid casting
    int minLambdaSerial = std::floor(numPerDex * log10(range.min()));
    int maxLambdaSerial = std::ceil(numPerDex * log10(range.max()));
    for (int k = minLambdaSerial; k <= maxLambdaSerial; ++k) lambdav.push_back(pow(10., k / numPerDex));

    // add the wavelengths mentioned in the configuration of the simulation
    for (double lambda : config->simulationWavelengths())
        if (range.contains(lambda)) lambdav.push_back(lambda);

    // add wavelength points around the threshold energies for all transitions
    for (const auto& pa_r : pa_resources)
    {
        for (double delta : {-2., -4. / 3., -2. / 3., 0., 2. / 3., 4. / 3., 2.})
        {
            double lambda = wavelengthToFromEnergy(pa_r.Eth + delta * pa_r.Es);
            if (range.contains(lambda)) lambdav.push_back(lambda);
        }
    }

    // add the fluorescence emission wavelengths
    for (const auto& fl_r : fl_resources)
    {
        double lambda = wavelengthToFromEnergy(fl_r.E);
        if (range.contains(lambda)) lambdav.push_back(lambda);
    }

    // add the outer wavelengths of our nonzero range, plus an extra just outside of that range,
    // so that there are always at least three points and thus two bins in the grid
    lambdav.push_back(nonZeroRange.min());
    lambdav.push_back(nonZeroRange.max());
    lambdav.push_back(nonZeroRange.max() * 1.000001);  // this wavelength point is never actually used

    // sort the wavelengths and remove duplicates
    NR::unique(lambdav);
    int numLambda = lambdav.size();

    // derive a wavelength grid that will be used for converting a wavelength to an index in the above array;
    // the grid points are shifted to the left of the actual sample points to approximate rounding
    _lambdav.resize(numLambda);
    _lambdav[0] = lambdav[0];
    for (int ell = 1; ell != numLambda; ++ell)
    {
        _lambdav[ell] = sqrt(lambdav[ell] * lambdav[ell - 1]);
    }

    // ---- extinction ---- //

    // calculate the extinction cross section at every wavelength; to guarantee that the cross section is zero
    // for wavelengths outside our range, leave the values for the three outer wavelength points at zero
    _sigmaextv.resize(numLambda);
    for (int ell = 1; ell < numLambda - 2; ++ell)
    {
        double lambda = lambdav[ell];
        double sigma = 0.;

        // bound electron scattering
        for (const auto& ion : _ions)
        {
            sigma += (_ray->sectionSca(lambda, ion.Z) + _com->sectionSca(lambda, ion.Z)) * ion.abund;
        }

        // photo-absorption and fluorescence
        for (const auto& pa_r : pa_resources)
        {
            const auto& ion = _ions[pa_r.ionIndex];
            double E = wavelengthToFromEnergy(lambda);
            sigma += pa_r.photoAbsorbThermalSection(E) * ion.abund;
        }

        // bound-bound absorption
        if (temperature() > 0.)
        {
            for (const auto& bt_res : bt_resources)
            {
                const auto& ion = _ions[bt_res.ionIndex];
                double un = _boundLevels[bt_res.upperIndex].pop;
                double ln = _boundLevels[bt_res.lowerIndex].pop;

                double Bul = bt_res.A * pow(bt_res.lam, 5.) / (2. * Constants::h() * Constants::c() * Constants::c());
                double Blu = Bul * bt_res.gul;

                double transrate = (ln * Blu - un * Bul) * ion.abund;
                if (transrate != 0.)
                {
                    constexpr double front = Constants::h() * Constants::c() / 4. / M_PI;

                    double stddev = bt_res.lam / Constants::c() * vtherm(ion.mass);
                    double gauss = gaussian(lambda, bt_res.lam, stddev);

                    sigma += front / bt_res.lam * transrate * gauss;
                }
            }
        }
        _sigmaextv[ell] = sigma;
    }

    // ---- scattering ---- //

    // make room for the scattering cross section and the cumulative fluorescence/scattering probabilities
    _sigmascav.resize(numLambda);
    _cumprobscavv.resize(numLambda, 0);

    // provide temporary array for the non-normalized fluorescence/scattering contributions (at the current wavelength)
    Array contribv(2 * _numIons + _numFluo);

    // calculate the above for every wavelength; as before, leave the values for the outer wavelength points at zero
    for (int ell = 1; ell < numLambda - 2; ++ell)
    {
        double lambda = lambdav[ell];
        double E = wavelengthToFromEnergy(lambda);

        // bound electron scattering
        for (int i = 0; i < _numIons; i++)
        {
            const auto& ion = _ions[i];

            contribv[i] = _ray->sectionSca(lambda, ion.Z) * _abundances[i];
            contribv[_numIons + i] = _com->sectionSca(lambda, ion.Z) * _abundances[i];
        }

        // fluorescence: iterate over both cross section and fluorescence parameter sets in sync
        for (int k = 0; k < _numFluo; k++)
        {
            const auto& flp_res = fl_resources[k];
            const auto& pap_res = pa_resources[flp_res.papIndex];
            double contribution = pap_res.photoAbsorbThermalSection(E) * _abundances[pap_res.ionIndex] * flp_res.omega;
            contribv[2 * _numIons + k] = contribution;
        }

        // determine the normalized cumulative probability distribution and the cross section
        _sigmascav[ell] = NR::cdf(_cumprobscavv[ell], contribv);
    }
}

////////////////////////////////////////////////////////////////////

XRayIonicGasMix::XRayIonicGasMix(SimulationItem* parent, string ionNames, BoundElectrons boundElectrons,
                                 vector<double> abundances, double temperature, bool setup)
{
    _ionNames = ionNames;
    _scatterBoundElectrons = boundElectrons;
    _abundances = abundances;
    _temperature = temperature;
    if (setup)
    {
        parent->addChild(this);
        this->setup();
    }
}

////////////////////////////////////////////////////////////////////

XRayIonicGasMix::~XRayIonicGasMix()
{
    delete _ray;
    delete _com;
}

////////////////////////////////////////////////////////////////////

int XRayIonicGasMix::indexForLambda(double lambda) const
{
    return NR::locateClip(_lambdav, lambda);
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::vtherm(double amu) const
{
    return sqrt(Constants::k() / Constants::amu() * temperature() / amu);
}

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType XRayIonicGasMix::materialType() const
{
    return MaterialMix::MaterialType::Gas;
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::hasPolarizedScattering() const
{
    return scatterBoundElectrons() == BoundElectrons::FreeWithPolarization;
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::hasExtraSpecificState() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::hasScatteringDispersion() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::hasLineEmission() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> XRayIonicGasMix::specificStateVariableInfo() const
{
    return vector<StateVariable>{StateVariable::numberDensity(), StateVariable::temperature()};
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::mass() const
{
    return Constants::Mproton();
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::sectionAbs(double lambda) const
{
    int index = indexForLambda(lambda);
    return _sigmaextv[index] - _sigmascav[index];
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::sectionSca(double lambda) const
{
    return _sigmascav[indexForLambda(lambda)];
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::sectionExt(double lambda) const
{
    return _sigmaextv[indexForLambda(lambda)];
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double number = state->numberDensity();
    return number > 0. ? sectionAbs(lambda) * number : 0.;
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::opacitySca(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double number = state->numberDensity();
    return number > 0. ? sectionSca(lambda) * number : 0.;
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double number = state->numberDensity();
    return number > 0. ? sectionExt(lambda) * number : 0.;
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMix::setScatteringInfoIfNeeded(PhotonPacket::ScatteringInfo* scatinfo, double lambda) const
{
    if (!scatinfo->valid)
    {
        scatinfo->valid = true;
        scatinfo->species = NR::locateClip(_cumprobscavv[indexForLambda(lambda)], random()->uniform());

        // for a fluorescence transition, determine the outgoing wavelength from the corresponding parameters
        if (scatinfo->species >= static_cast<int>(2 * _numIons))
        {
            int k = scatinfo->species - 2 * _numIons;
            const auto& flp = _fluorescence[k];
            if (temperature() > 0.)
            {
                const auto& ion = _ions[flp.ionIndex];
                scatinfo->velocity = ion.vth * random()->maxwell();
            }

            if (!flp.W)
            {
                // for a zero-width line, simply copy the central wavelength
                scatinfo->lambda = flp.E;  // already converted to wavelength
            }
            else
            {
                // otherwise sample a wavelength from the Lorentz line shape in energy space;
                // the tails of the Lorentz distribition are very long, occasionaly resulting in negative energies;
                // therefore we loop until the sampled wavelength is meaningful
                while (true)
                {
                    // fix this!! use a biasing technique here?
                    scatinfo->lambda = wavelengthToFromEnergy(flp.E + flp.W * random()->lorentz());
                    if (nonZeroRange.contains(scatinfo->lambda)) break;
                }
            }
        }
        else
        {
            int i = scatinfo->species % _numIons;
            const auto& ion = _ions[i];
            scatinfo->velocity = ion.vth * random()->maxwell();
        }
    }
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMix::peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs,
                                        Direction bfky, const MaterialState* /*state*/, const PhotonPacket* pp) const
{
    // draw a random scattering channel and atom velocity, unless a previous peel-off stored this already
    auto scatinfo = const_cast<PhotonPacket*>(pp)->getScatteringInfo();
    setScatteringInfoIfNeeded(scatinfo, lambda);

    // if we have dispersion, for electron scattering, adjust the incoming wavelength to the electron rest frame
    if (temperature() > 0. && scatinfo->species < static_cast<int>(2 * _numIons))
        lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);

    // Rayleigh scattering in electron rest frame; no support for polarization
    if (scatinfo->species < static_cast<int>(_numIons))
    {
        _ray->peeloffScattering(I, lambda, scatinfo->species + 1, pp->direction(), bfkobs);
    }

    // Compton scattering in electron rest frame; with support for polarization if enabled
    else if (scatinfo->species < static_cast<int>(2 * _numIons))
    {
        _com->peeloffScattering(I, Q, U, V, lambda, scatinfo->species - _numIons + 1, pp->direction(), bfkobs, bfky,
                                pp);
    }

    // fluorescence
    else
    {
        // unpolarized isotropic emission; the bias weight is trivially 1 and there is no contribution to Q, U, V
        I = 1.;

        // update the photon packet wavelength to the (possibly sampled) wavelength of this fluorescence transition
        lambda = scatinfo->lambda;
    }

    // if we have dispersion, Doppler-shift the outgoing wavelength from the electron rest frame
    if (temperature() > 0.) lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfkobs, scatinfo->velocity);
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMix::performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const
{
    // draw a random fluorescence channel and atom velocity, unless a previous peel-off stored this already
    auto scatinfo = pp->getScatteringInfo();
    setScatteringInfoIfNeeded(scatinfo, lambda);

    // if we have dispersion, for electron scattering, adjust the incoming wavelength to the electron rest frame
    if (temperature() > 0. && scatinfo->species < static_cast<int>(2 * _numIons))
        lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);

    // room for the outgoing direction
    Direction bfknew;

    // Rayleigh scattering, no support for polarization: determine the new propagation direction
    if (scatinfo->species < static_cast<int>(_numIons))
    {
        bfknew = _ray->performScattering(lambda, scatinfo->species + 1, pp->direction());
    }

    // Compton scattering, with support for polarization if enabled:
    // determine the new propagation direction and wavelength, and if polarized, update the stokes vector
    else if (scatinfo->species < static_cast<int>(2 * _numIons))
    {
        bfknew = _com->performScattering(lambda, scatinfo->species - _numIons + 1, pp->direction(), pp);
    }

    // fluorescence, always unpolarized and isotropic
    else
    {
        // update the photon packet wavelength to the (possibly sampled) wavelength of this fluorescence transition
        lambda = scatinfo->lambda;

        // draw a random, isotropic outgoing direction
        bfknew = random()->direction();

        // clear the stokes vector (only relevant if polarization support is enabled)
        pp->setUnpolarized();
    }

    // if we have dispersion, Doppler-shift the outgoing wavelength from the electron rest frame
    if (temperature() > 0.) lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfknew, scatinfo->velocity);

    // execute the scattering event in the photon packet
    pp->scatter(bfknew, state->bulkVelocity(), lambda);
}

////////////////////////////////////////////////////////////////////

Array XRayIonicGasMix::lineEmissionCenters() const
{
    Array centers(_numBTransitions);
    for (int b = 0; b != _numBTransitions; ++b)
    {
        const auto& bbp = _boundTransitions[b];
        centers[b] = bbp.lam;
    }
    return centers;
}

////////////////////////////////////////////////////////////////////

Array XRayIonicGasMix::lineEmissionMasses() const
{
    Array mass(_numBTransitions);
    for (int b = 0; b != _numBTransitions; ++b)
    {
        const auto& bbp = _boundTransitions[b];
        const auto& ion = _ions[bbp.ionIndex];
        mass[b] = ion.mass;
    }
    return mass;
}

////////////////////////////////////////////////////////////////////

Array XRayIonicGasMix::lineEmissionSpectrum(const MaterialState* state, const Array& /*Jv*/) const
{
    Array luminosities(_numBTransitions);
    double front = Constants::h() * Constants::c() * state->volume();
    for (int b = 0; b < _numBTransitions; b++)
    {
        const auto& bt = _boundTransitions[b];
        const auto& bl = _boundLevels[bt.upperIndex];
        luminosities[b] = front / bt.lam * bt.A * bl.pop * state->numberDensity();
    }
    return luminosities;
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::indicativeTemperature(const MaterialState* /*state*/, const Array& /*Jv*/) const
{
    return temperature();
}

////////////////////////////////////////////////////////////////////
