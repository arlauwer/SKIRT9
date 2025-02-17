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

////////////////////////////////////////////////////////////////////

struct FluorescenceParams
{
    FluorescenceParams() {}

    FluorescenceParams(const Array& a) : omega(a[0]), E(a[1]), W(a[2]) {}

    // short Z;       // atomic number
    // short n;       // principal quantum number of the shell
    // short l;       // orbital quantum number of the subshell
    double omega;  // fluorescence yield (1)
    double E;      // (central) energy of the emitted photon (eV)
    double W;      // FWHM of the Lorentz shape for the emitted photon (eV), or zero
};

struct PhotoAbsorbParams
{
    PhotoAbsorbParams() {}

    PhotoAbsorbParams(const Array& a)
        : Z(a[0]), n(a[1]), l(a[2]), Eth(a[3]), Emax(5e5), E0(a[4]), sigma0(a[5]), ya(a[6]), P(a[7]), yw(a[8]), y0(0.),
          y1(0.)
    {}

    vector<FluorescenceParams> fluorescenceParamsList;

    short Z;        // atomic number
    short n;        // principal quantum number of the shell
    short l;        // orbital quantum number of the subshell
    double Eth;     // subshell ionization threshold energy (eV)
    double Emax;    // maximum energy for validity of the formula (eV)
    double E0;      // fit parameter (eV)
    double sigma0;  // fit parameter (Mb = 10^-22 m^2)
    double ya;      // fit parameter (1)
    double P;       // fit parameter (1)
    double yw;      // fit parameter (1)
    double y0;      // fit parameter (1)
    double y1;      // fit parameter (1)
};

struct IonParams
{
    IonParams(const Array& a) : Z(a[0]), N(a[1]), mass(a[2]) {}

    short Z;
    short N;
    double mass;

    PhotoAbsorbParams photoAbsorbParams;
};

namespace
{
    // ---- common helper functions ----

    // convert photon energy in eV to and from wavelength in m (same conversion in both directions)
    constexpr double wavelengthToFromEnergy(double x)
    {
        constexpr double front = Constants::h() * Constants::c() / Constants::Qelectron();
        return front / x;
    }

    // return thermal velocity for given gas temperature (in K) and particle mass (in amu)
    double vtherm(double T, double amu)
    {
        return sqrt(Constants::k() / Constants::amu() * T / amu);
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

    // ---- hardcoded configuration constants ----

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
    template<class S, int N> vector<S*> loadStruct(const SimulationItem* item, string filename, string description)
    {
        vector<S*> result;
        TextInFile infile(item, filename, description, true);
        for (int i = 0; i != N; ++i) infile.addColumn(string());
        Array row;
        while (infile.readRow(row)) result.emplace_back(new S(row));
        return result;
    }

    // ---- photo-absorption cross section ----

    // return photo-absorption cross section in m2 for given energy in eV and cross section parameters,
    // without taking into account thermal dispersion
    double crossSection(double E, const PhotoAbsorbParams& p)
    {
        if (E < p.Eth || E >= p.Emax) return 0.;

        double x = E / p.E0 - p.y0;
        double y = std::sqrt(x * x + p.y1 * p.y1);
        double xm1 = x - 1.;
        double Q = 5.5 + p.l - 0.5 * p.P;
        double F = (xm1 * xm1 + p.yw * p.yw) * std::pow(y, -Q) * std::pow(1. + std::sqrt(y / p.ya), -p.P);
        return 1e-22 * p.sigma0 * F;  // from Mb to m2
    }

    // return photo-absorption cross section in m2 for given energy in eV and cross section parameters,
    // approximating thermal dispersion by replacing the steep threshold transition by a sigmoid
    // error function with given parameters (dispersion and maximum value)
    double crossSection(double E, double Es, double sigmamax, const PhotoAbsorbParams& p)
    {
        if (E <= p.Eth - 2. * Es) return 0.;
        if (E >= p.Eth + 2. * Es) return crossSection(E, p);
        return sigmamax * (0.5 + 0.5 * std::erf((E - p.Eth) / Es));
    }

    // ---- bound-electron scattering resources ----

    // load data from resource file with N columns into a vector of N arrays, and return that vector;
    // each of the arrays is resized to remove trailing NaN values, if applicable
    vector<Array> loadColumns(int N, const SimulationItem* item, string filename, string description)
    {
        TextInFile infile(item, filename, description, true);
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
}

////////////////////////////////////////////////////////////////////

// ---- base class for scattering helpers ----

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

// ---- no scattering helper ----

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

// ---- free-electron Compton scattering helper ----

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

// ---- free-electron Compton with polarization scattering helper ----

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

// ---- bound-electron Compton scattering helper ----

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
        BoundComptonHelper(SimulationItem* item, int numAtoms)
        {
            // load bound Compton cross sections
            _CSv = loadColumns(numAtoms + 1, item, "XRay_CS.txt", "bound Compton data");
            _CSv[0] *= keVtoScaledEnergy;                            // convert from keV to 1
            for (size_t Z = 1; Z <= numAtoms; ++Z) _CSv[Z] *= 1e-4;  // convert from cm2 to m2

            // load incoherent scattering functions
            _SFv = loadColumns(numAtoms + 1, item, "XRay_SF.txt", "bound Compton data");

            // load pdfs for projected momentum of target electron
            _CPv = loadColumns(numAtoms + 1, item, "XRay_CP.txt", "bound Compton data");
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
            for (size_t Z = 2; Z <= numAtoms; ++Z)
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

// ---- smooth Rayleigh scattering helper ----

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
        SmoothRayleighHelper(SimulationItem* item, int numAtoms)
        {
            // load smooth Rayleigh cross sections
            _RSSv = loadColumns(numAtoms + 1, item, "XRay_RSS.txt", "smooth Rayleigh data");
            _RSSv[0] *= keVtoScaledEnergy;                         // convert from keV to 1
            for (int Z = 1; Z <= numAtoms; ++Z) _RSSv[Z] *= 1e-4;  // convert from cm2 to m2

            // load atomic form factors
            _FFv = loadColumns(numAtoms + 1, item, "XRay_FF.txt", "smooth Rayleigh data");

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

// ---- anomalous Rayleigh scattering helper ----

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
        AnomalousRayleighHelper(SimulationItem* item, int numAtoms)
        {
            // load anomalous Rayleigh cross sections, atomic form factors and anomalous scattering functions
            _RSAv = loadColumns(2 * numAtoms + 2, item, "XRay_RSA.txt", "anomalous Rayleigh data");
            _FFv = loadColumns(numAtoms + 1, item, "XRay_FF.txt", "anomalous Rayleigh data");
            _F1v = loadColumns(2 * numAtoms + 2, item, "XRay_F1.txt", "anomalous Rayleigh data");
            _F2v = loadColumns(2 * numAtoms + 2, item, "XRay_F2.txt", "anomalous Rayleigh data");

            // convert units
            for (size_t Z = 1; Z <= numAtoms; ++Z)
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

    // ---- load resources ----

    // load all ions WIP
    _numIons = 1;
    _numAtoms = 1;
    _numFluos = 1;
    _ionParams.push_back(new IonParams({1, 1, 1.}));
    _ionParams[0]->photoAbsorbParams = PhotoAbsorbParams({1, 1, 0, 13.60, 0.4298, 54750, 32.88, 2.963, 0});
    _ionParams[0]->photoAbsorbParams.fluorescenceParamsList.push_back(FluorescenceParams({1, 5, 0}));
    ///////// WIP

    _fluorescenceMass.resize(_numFluos);
    _fluorescenceParams.resize(_numFluos);
    int f = 0;
    for (int i = 0; i < _numIons; ++i)
    {
        const auto ion = _ionParams[i];
        const auto& pap = ion->photoAbsorbParams;
        for (const auto& flp : pap.fluorescenceParamsList)
        {
            _fluorescenceMass[f] = ion->mass;
            _fluorescenceParams[f] = &flp;
            f++;
        }
    }

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
            _ray = new SmoothRayleighHelper(this, _numAtoms);
            _com = new BoundComptonHelper(this, _numAtoms);
            break;
        case BoundElectrons::Exact:
            _ray = new AnomalousRayleighHelper(this, _numAtoms);
            _com = new BoundComptonHelper(this, _numAtoms);
            break;
    }
}

////////////////////////////////////////////////////////////////////

XRayIonicGasMix::~XRayIonicGasMix()
{
    delete _ray;
    delete _com;
}

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType XRayIonicGasMix::materialType() const
{
    return MaterialMix::MaterialType::Gas;
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::hasPolarizedScattering() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::hasExtraSpecificState() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::hasScatteringDispersion() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> XRayIonicGasMix::parameterInfo() const
{
    vector<SnapshotParameter> descriptors;
    descriptors.push_back(SnapshotParameter::temperature());

    for (int f = 0; f != _numIons; ++f)
    {
        descriptors.push_back(SnapshotParameter::custom("ion"));
    }
    return descriptors;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> XRayIonicGasMix::specificStateVariableInfo() const
{
    // add hard-coded ions here for now. Later we should try to use some sort of useColumns="auto" to automatically
    // determine the columns
    vector<StateVariable> vars = {StateVariable::numberDensity(), StateVariable::temperature()};
    for (int i = 0; i < _numIons; ++i)
    {
        vars.push_back(StateVariable::custom(i, "ion" + std::to_string(i), "ionization fraction"));
    }
    return vars;
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMix::initializeSpecificState(MaterialState* state, double metallicity, double temperature,
                                              const Array& params) const
{
    if (state->numberDensity() > 0.)
    {
        double defaultTemperature = 1;
        state->setTemperature(temperature >= 0. ? temperature : defaultTemperature);

        for (int i = 0; i < _numIons; ++i)
        {
            state->setCustom(i, 1.); // hard code the ionization fraction to 1 for now
        }
    }
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::mass() const
{
    return Constants::Mproton();
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::sectionAbs(double /*lambda*/) const
{
    return 0.0;
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::sectionSca(double /*lambda*/) const
{
    return 0.0;
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::sectionExt(double /*lambda*/) const
{
    return 0.0;
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const
{
    return opacityExt(lambda, state, pp) - opacitySca(lambda, state, pp);
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::opacitySca(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    return sigmaSca(lambda, state).sum();
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double sigma = 0.;
    double E = wavelengthToFromEnergy(lambda);

    // bound electron scattering
    for (int i = 0; i < _numIons; ++i)
    {
        double abund = state->custom(i);
        int Z = _ionParams[i]->Z;

        sigma += (_ray->sectionSca(lambda, Z) + _com->sectionSca(lambda, Z)) * abund;
    }

    // photo-absorption and fluorescence
    for (int i = 0; i < _numIons; ++i)
    {
        const auto ion = _ionParams[i];
        const auto pap = ion->photoAbsorbParams;
        double mass = ion->mass;
        double abund = state->custom(i);

        // use thermal sigmoid
        double Es = pap.Eth * vtherm(state->temperature(), mass) / Constants::c();
        double sigmamax = crossSection(pap.Eth + 2. * Es, pap);

        sigma += crossSection(E, Es, sigmamax, pap) * abund;
    }

    double number = state->numberDensity();
    return number > 0. ? sigma * number : 0.;
}

////////////////////////////////////////////////////////////////////

Array XRayIonicGasMix::sigmaSca(double lambda, const MaterialState* state) const
{
    Array sigma(2 * _numIons + _numFluos);
    double E = wavelengthToFromEnergy(lambda);

    // bound electron scattering
    for (int i = 0; i < _numIons; ++i)
    {
        double abund = state->custom(i);
        int Z = _ionParams[i]->Z;

        sigma[i] = _ray->sectionSca(lambda, Z) * abund;
        sigma[_numIons + i] = _com->sectionSca(lambda, Z) * abund;
    }

    // fluorescence
    int f = 0;
    for (int i = 0; i < _numIons; ++i)
    {
        const auto ion = _ionParams[i];
        const auto pap = ion->photoAbsorbParams;
        double mass = ion->mass;
        double abund = state->custom(i);

        // use thermal sigmoid
        double Es = pap.Eth * vtherm(state->temperature(), mass) / Constants::c();
        double sigmamax = crossSection(pap.Eth + 2. * Es, pap);

        for (const auto& flp : pap.fluorescenceParamsList)
        {
            sigma[2 * _numIons + f] = crossSection(E, Es, sigmamax, pap) * abund * flp.omega;
            f++;
        }
    }

    return sigma;
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMix::setScatteringInfoIfNeeded(PhotonPacket::ScatteringInfo* scatinfo, double lambda,
                                                const MaterialState* state) const
{
    if (!scatinfo->valid)
    {
        scatinfo->valid = true;

        // determine the scattering channel
        Array sigma = sigmaSca(lambda, state);
        Array cumsigma;
        NR::cdf(cumsigma, sigma);
        scatinfo->species = NR::locateClip(cumsigma, random()->uniform());

        // determine thermal velocity
        if (state->temperature() > 0)
        {
            double v;
            if (scatinfo->species < 2 * _numIons)  // bound electron scattering
                v = vtherm(state->temperature(), _ionParams[scatinfo->species % _numIons]->mass);
            else  // fluorescence
                v = vtherm(state->temperature(), _fluorescenceMass[scatinfo->species - 2 * _numIons]);
            scatinfo->velocity = v * random()->maxwell();
        }

        // for a fluorescence transition, determine the outgoing wavelength from the corresponding parameters
        if (scatinfo->species >= 2 * _numIons)
        {
            int f = scatinfo->species - 2 * _numIons;
            const auto* flp = _fluorescenceParams[f];

            if (!flp->W)
            {
                // if the line shape has zero width, convert the central enery to the fixed wavelength
                scatinfo->lambda = wavelengthToFromEnergy(flp->E);
            }
            else
            {
                while (true)
                {
                    scatinfo->lambda = wavelengthToFromEnergy(flp->E + flp->W / 2. * random()->lorentz());
                    if (nonZeroRange.contains(scatinfo->lambda)) break;
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMix::peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs,
                                        Direction bfky, const MaterialState* state, const PhotonPacket* pp) const
{
    // draw a random scattering channel and atom velocity, unless a previous peel-off stored this already
    auto scatinfo = const_cast<PhotonPacket*>(pp)->getScatteringInfo();
    setScatteringInfoIfNeeded(scatinfo, lambda, state);

    // if we have dispersion, for electron scattering, adjust the incoming wavelength to the electron rest frame
    if (state->temperature() > 0. && scatinfo->species < 2 * _numIons)
        lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);

    // Rayleigh scattering in electron rest frame; no support for polarization
    if (scatinfo->species < _numIons)
    {
        const auto ion = _ionParams[scatinfo->species];
        _ray->peeloffScattering(I, lambda, ion->Z, pp->direction(), bfkobs);
    }

    // Compton scattering in electron rest frame; with support for polarization if enabled
    else if (scatinfo->species < 2 * _numIons)
    {
        const auto ion = _ionParams[scatinfo->species - _numIons];
        _com->peeloffScattering(I, Q, U, V, lambda, ion->Z, pp->direction(), bfkobs, bfky, pp);
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
    if (state->temperature() > 0.) lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfkobs, scatinfo->velocity);
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMix::performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const
{
    // draw a random fluorescence channel and atom velocity, unless a previous peel-off stored this already
    auto scatinfo = pp->getScatteringInfo();
    setScatteringInfoIfNeeded(scatinfo, lambda, state);

    // if we have dispersion, for electron scattering, adjust the incoming wavelength to the electron rest frame
    if (state->temperature() > 0. && scatinfo->species < 2 * _numIons)
        lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);

    // room for the outgoing direction
    Direction bfknew;

    // Rayleigh scattering, no support for polarization: determine the new propagation direction
    if (scatinfo->species < _numIons)
    {
        const auto ion = _ionParams[scatinfo->species];
        bfknew = _ray->performScattering(lambda, ion->Z, pp->direction());
    }

    // Compton scattering, with support for polarization if enabled:
    // determine the new propagation direction and wavelength, and if polarized, update the stokes vector
    else if (scatinfo->species < 2 * _numIons)
    {
        const auto ion = _ionParams[scatinfo->species - _numIons];
        bfknew = _com->performScattering(lambda, ion->Z, pp->direction(), pp);
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
    if (state->temperature()) lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfknew, scatinfo->velocity);

    // execute the scattering event in the photon packet
    pp->scatter(bfknew, state->bulkVelocity(), lambda);
}
