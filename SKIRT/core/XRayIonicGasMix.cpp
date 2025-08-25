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
#include "ListWavelengthGrid.hpp"
#include "MaterialState.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "Range.hpp"
#include "StringUtils.hpp"
#include <map>
#include <regex>

////////////////////////////////////////////////////////////////////

namespace
{
    // ---- common helper functions ---- //

    static constexpr double defaultTemperature = 1.;
    static constexpr double defaultMetallicity = 0.;
    // The amount of bins in the radiation field (these should match the Cloudy table)
    static constexpr int radBins = 2;

    static constexpr int numAtoms = 30;  // H to Zn
    // masses of the elements in amu up to Z=30
    const vector<double> atomicMasses = {1.0079,  4.0026, 6.941,   9.01218, 10.81,   12.011,  14.0067, 15.9994,
                                         18.9984, 20.179, 22.9898, 24.305,  26.9815, 28.0855, 30.9738, 32.06,
                                         35.453,  39.948, 39.0983, 40.08,   44.9559, 47.9,    50.9415, 51.996,
                                         54.938,  55.847, 58.9332, 58.7,    63.546,  65.38};

    // we use a struct to organize Z, N instead of using a function which would require a sqrt, etc.
    struct Ion
    {
        short Z;  // atomic number
        short N;  // number of electrons
    };

    static const vector<Ion> initIons()
    {
        vector<Ion> ions;
        for (short Z = 1; Z <= numAtoms; ++Z)
        {
            for (short N = 1; N <= Z; ++N)
            {
                Ion ion;
                ion.Z = Z;
                ion.N = N;
                ions.push_back(ion);
            }
        }
        return ions;
    }

    // organized way to order the ions by Z and N

    static const vector<Ion> ions = initIons();  // indexed on ions
    static const int numIons = ions.size();      // total number of ions (i.e. abundances.size())
    // WIP THE NUMIONS WILL PROBABLY HAVE TO BE CHANGED? IF THE ENTIRE CLOUDY TABLE DOESN'T HAVE EG. Zn+24
    // THEN JUST DON'T STORE IT AND numIons < 465

    double vtherm(double T, double amu)
    {
        return sqrt(Constants::k() / Constants::amu() * T / amu);
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

    // ---- bound-electron scattering resources ---- //

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
        vector<Array> _CSv;  // 0: E (keV->1); 1-numAtoms: bound Compton cross sections (cm2->m2)
        vector<Array> _SFv;  // 0: q (1); 1-numAtoms: incoherent scattering functions (1)
        vector<Array> _CPv;  // 0: E (keV->1); 1-numAtoms: pdf for target electron momentum (1)
        vector<Array> _IBv;  // 0: E (keV->1) ionisation energy of the outer subshell electrons

        // precalculated cumulative distributions for target electron momentum
        Range _cumRange;
        vector<Array> _cumCPv;  // 0: E axis; 1-numAtoms: cumulative pdf for target electron momentum

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

// ---- smooth Rayleigh scattering helper ---- //

namespace
{
    // this helper implements smooth Rayleigh scattering;
    // below the energy limit of the tabulated data, use Thomson scattering instead
    class SmoothRayleighHelper : public XRayIonicGasMix::ScatteringHelper
    {
    private:
        vector<Array> _RSSv;  // 0: E (keV->1); 1-numAtoms: smooth Rayleigh cross sections (cm2->m2)
        vector<Array> _FFv;   // 0: q (1); 1-numAtoms: atomic form factors (1)
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

// ---- anomalous Rayleigh scattering helper ---- //

namespace
{
    // this helper implements anomalous Rayleigh scattering
    // below the energy limit of the tabulated data, use Thomson scattering instead
    class AnomalousRayleighHelper : public XRayIonicGasMix::ScatteringHelper
    {
    private:
        vector<Array> _RSAv;  // 2*Z: E (keV->1); 2*Z+1: anomalous Rayleigh cross sections (cm2->m2)
        vector<Array> _FFv;   // 0: q (1); 1-numAtoms: atomic form factors (1)
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

    // ---- load optical properties ---- //
    _abundanceTable.open(this, "abund.stab", "ion(1),n(1/m3),Z(1),bin0(W/m3),bin1(W/m3)", "abund(1/m3)", true, false);
    _temperatureTable.open(this, "temp.stab", "n(1/m3),Z(1),bin0(W/m3),bin1(W/m3)", "temp(K)", true, false);

    _opacityTable.open(this, "opt.stab", "lam(m),n(1/m3),Z(1),bin0(W/m3),bin1(W/m3)", "opac(1/m)", true, false);
    _emissivityTable.open(this, "opt.stab", "lam(m),n(1/m3),Z(1),bin0(W/m3),bin1(W/m3)", "emis(W/m3)", true, false);

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

    // lambda wavelength grid
    _range = config()->simulationWavelengthRange();
    _range.intersect(nonZeroRange);

    // slice the radiation field to the configured wavelength range
    Array rad = config()->radiationFieldWLG()->lambdav();
    if (rad.size() != radBins)
    {
        throw FATALERROR(
            "The radiation field wavelength grid must have the same number of bins as the precomputed table");
    }

    // obtain the wavelengths from the absorption table
    _opacityTable.axisArray<0, double>(_lambda);
    _numLambda = _lambda.size();

    if (_lambda[0] > _lambda[_numLambda - 1])
    {
        throw FATALERROR("The tabulated wavelength grid must be ascending");
    }

    // ---- continuum emissivity ---- //
    vector<double> lambdaCv(std::begin(_lambda), std::end(_lambda));
    if (!_emissionGrid) _emissionGrid = new ListWavelengthGrid(this, lambdaCv, 0., true, true);

    // add 2 outer points to match the _emissionGrid
    // PREFERABLY DO THIS IN THE CLOUDY TABLE BEFOREHAND
    // Array temp = _emissivityC;
    // _emissivityC.resize(_numLambda + 2);
    // _emissivityC[0] = 0.;
    // _emissivityC[std::slice(1, _numLambda, 1)] = temp;
    // _emissivityC[_numLambda + 1] = 0.;
}

////////////////////////////////////////////////////////////////////

XRayIonicGasMix::~XRayIonicGasMix()
{
    delete _ray;
    delete _com;
    delete _emissionGrid;
}

////////////////////////////////////////////////////////////////////

int XRayIonicGasMix::indexForLambda(double lambda) const
{
    return NR::locateClip(_lambda, lambda);
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

MaterialMix::DynamicStateType XRayIonicGasMix::hasDynamicMediumState() const
{
    return DynamicStateType::Primary;
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::hasScatteringDispersion() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::hasContinuumEmission() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::hasLineEmission() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> XRayIonicGasMix::specificStateVariableInfo() const
{
    vector<StateVariable> result{StateVariable::numberDensity(), StateVariable::temperature(),
                                 StateVariable::metallicity()};

    // next available custom variable index
    int index = 0;

    // State variables
    // abundances:      numIons
    // vtherm:          numAtoms
    // sigmaabsC:       numLambda
    // sigmascaC:       numLambda
    // cumprobscaC:     numLambda x 2*numIons+1 (+1 for cumulative)
    // emissivityC:     numLambda + 2
    // LINES WIP

    const_cast<XRayIonicGasMix*>(this)->_indexAbundances = index;
    for (int Z = 1; Z <= numAtoms; ++Z)
        for (int N = 0; N != Z; ++N)
            result.push_back(StateVariable::custom(index++, "abundance", "numbervolumedensity"));

    const_cast<XRayIonicGasMix*>(this)->_indexThermalVelocity = index;
    result.push_back(StateVariable::custom(index++, "thermal velocity", "velocity"));

    const_cast<XRayIonicGasMix*>(this)->_indexSigmaAbs = index;
    for (int l = 0; l < _numLambda; l++) result.push_back(StateVariable::custom(index++, "absorption opacity", "area"));

    const_cast<XRayIonicGasMix*>(this)->_indexSigmaSca = index;
    for (int l = 0; l < _numLambda; l++) result.push_back(StateVariable::custom(index++, "scattering opacity", "area"));

    const_cast<XRayIonicGasMix*>(this)->_indexSigmaScaCum = index;
    for (int l = 0; l < _numLambda; l++)
        for (int i = 0; i < 2 * numIons + 1; ++i)
            result.push_back(StateVariable::custom(index++, "cumulative scattering probability", "1"));

    const_cast<XRayIonicGasMix*>(this)->_indexEmissivity = index;
    for (int l = 0; l < _numLambda + 2; l++)
        result.push_back(StateVariable::custom(index++, "volume emissivity", "powervolumedensity"));

    return result;
}

////////////////////////////////////////////////////////////////////

#define setAbundance(ion, value) setCustom(_indexAbundances + (ion), (value))
#define getAbundance(ion) custom(_indexAbundances + (ion))
#define setVTherm(Z, value) setCustom(_indexThermalVelocity + (Z), (value))
#define getVTherm(Z) custom(_indexThermalVelocity + (Z))
#define setSigmaAbs(l, value) setCustom(_indexSigmaAbs + (l), (value))
#define getSigmaAbs(l) custom(_indexSigmaAbs + (l))
#define setSigmaSca(l, value) setCustom(_indexSigmaSca + (l), (value))
#define getSigmaSca(l) custom(_indexSigmaSca + (l))
#define setSigmaScaCum(l, index, value) setCustom(_indexSigmaScaCum + (l) * 2 * numAtoms + (index), (value))
#define getSigmaScaCum(l, index) custom(_indexSigmaScaCum + (l) * 2 * numAtoms + (index))
#define setEmissivity(l, value) setCustom(_indexEmissivity + (l), (value))
#define getEmissivity(l) custom(_indexEmissivity + (l))

////////////////////////////////////////////////////////////////////

void XRayIonicGasMix::initializeSpecificState(MaterialState* state, double Z, double /*temperature*/,
                                              const Array& /*params*/) const
{
    double J1 = 0.;  // determine this using lum and grid extent?
    double J2 = 0.;  // determine this using lum and grid extent?
    double n = state->numberDensity();
    Z = Z >= 0. ? Z : defaultMetallicity;

    double T = _temperatureTable(n, Z, J1, J2);

    state->setTemperature(T);
    state->setMetallicity(Z);

    for (int i = 0; i < numIons; i++)
    {
        double abund = _abundanceTable((double)i, n, Z, J1, J2);
        state->setAbundance(i, abund);
    }

    for (int Z = 0; Z < numAtoms; Z++)
    {
        double v = vtherm(T, atomicMasses[Z]);
        state->setVTherm(Z, v);
    }

    for (int l = 0; l < _numLambda; l++)
    {
        double lambda = _lambda[l];

        double abs = _opacityTable((double)l, n, Z, J1, J2);
        double emi = _emissivityTable((double)l, n, Z, J1, J2);

        state->setSigmaAbs(l, abs);
        state->setEmissivity(l, emi);

        // Scattering
        // provide temporary array for the non-normalized fluorescence/scattering contributions (at the current wavelength)
        Array sigmaScaFractions(2 * numIons);
        Array sigmaScaCum;

        // bound electron scattering
        for (int i = 0; i < numIons; i++)
        {
            double abund = state->getAbundance(i);

            int Z = ions[i].Z;  // TEMP THIS NEEDS TO BE FIXED LIKE IN CLOUDY WHERE FREE ELECTRON FRACTION IS USED

            sigmaScaFractions[i] = _ray->sectionSca(lambda, Z) * abund;
            sigmaScaFractions[numIons + i] = _com->sectionSca(lambda, Z) * abund;
        }

        // determine the normalized cumulative probability distribution and the cross section
        double sigmaSca = NR::cdf(sigmaScaCum, sigmaScaFractions);

        state->setSigmaSca(l, sigmaSca);
        for (int i = 0; i < 2 * numIons + 1; i++)
        {
            state->setSigmaScaCum(l, i, sigmaScaCum[i]);
        }
    }
}

////////////////////////////////////////////////////////////////////

UpdateStatus XRayIonicGasMix::updateSpecificState(MaterialState* state, const Array& Jv) const
{
    UpdateStatus status;

    Array radWidth = config()->radiationFieldWLG()->dlambdav();

    // Lookup table values
    Array J = Jv * radWidth * 4. * M_PI;  // W/m2/m/sr -> W/m2 (integrated mean intensity)
    double J1 = J[0];
    double J2 = J[1];
    double n = state->numberDensity();
    double Z = state->metallicity();

    Array prevAbund(numIons);

    for (int i = 0; i < numIons; i++)
    {
        prevAbund[i] = state->getAbundance(i);
        double abund = _abundanceTable((double)i, n, Z, J1, J2);
        state->setAbundance(i, abund);
    }

    double temp = _temperatureTable(n, Z, J1, J2);
    state->setTemperature(temp);

    for (int i = 0; i < numAtoms; i++)
    {
        double v = vtherm(temp, atomicMasses[i]);
        state->setVTherm(i, v);
    }

    for (int l = 0; l < _numLambda; l++)
    {
        double lambda = _lambda[l];

        double abs = _opacityTable(lambda, n, Z, J1, J2);
        double emi = _emissivityTable(lambda, n, Z, J1, J2);

        state->setSigmaAbs(l, abs);
        state->setEmissivity(l, emi);

        // Scattering
        // provide temporary array for the non-normalized fluorescence/scattering contributions (at the current wavelength)
        Array sigmaScaFractions(2 * numIons);
        Array sigmaScaCum;

        // bound electron scattering
        for (int i = 0; i < numIons; i++)
        {
            double abund = state->getAbundance(i);

            int Z = ions[i].Z;  // TEMP THIS NEEDS TO BE FIXED LIKE IN CLOUDY WHERE FREE ELECTRON FRACTION IS USED

            sigmaScaFractions[i] = _ray->sectionSca(lambda, Z) * abund;
            sigmaScaFractions[numIons + i] = _com->sectionSca(lambda, Z) * abund;
        }

        // determine the normalized cumulative probability distribution and the cross section
        double sigmaSca = NR::cdf(sigmaScaCum, sigmaScaFractions);

        state->setSigmaSca(l, sigmaSca);
        for (int i = 0; i < 2 * numIons + 1; i++)
        {
            state->setSigmaScaCum(l, i, sigmaScaCum[i]);
        }
    }

    // calculate the amount of cells that were clipped (i.e. not interpolated in the lookup table)
    double conv = 0.;
    for (int i = 0; i < numIons; i++)
    {
        double diff = (prevAbund[i] - state->getAbundance(i));
        conv += diff * diff;
    }
    // somehow handle the zero-abundances

    // status.updateConverged();
    status.updateNotConverged();
    return status;
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::isSpecificStateConverged(int numCells, int numUpdated, int numNotConverged,
                                               MaterialState* /*currentAggregate*/,
                                               MaterialState* /*previousAggregate*/) const
{
    return numCells == numUpdated && numNotConverged == 0;
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::mass() const
{
    return Constants::Mproton();
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::sectionAbs(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::sectionSca(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::sectionExt(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    int l = indexForLambda(lambda);
    return state->getSigmaAbs(l);
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::opacitySca(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    int l = indexForLambda(lambda);
    return state->getSigmaSca(l);
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    return opacityAbs(lambda, state, nullptr) + opacitySca(lambda, state, nullptr);
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMix::setScatteringInfoIfNeeded(PhotonPacket::ScatteringInfo* scatinfo, double lambda,
                                                const MaterialState* state) const
{
    if (!scatinfo->valid)
    {
        scatinfo->valid = true;
        int lam = NR::locateClip(_lambda, lambda);  // this should (almost) never have to clip
        // scattering can only happen if opacity is non-zero, so lambda should be in range of _lambdaC
        // maybe some Doppler shift but a simple clip should be sufficient
        Array sigmaScaCum(2 * numIons + 1);
        for (int i = 0; i < 2 * numIons + 1; i++) sigmaScaCum[i] = state->getSigmaScaCum(lam, i);

        scatinfo->species = NR::locateClip(sigmaScaCum, random()->uniform());

        int i = scatinfo->species % numIons;
        int Z = ions[i].Z;

        scatinfo->velocity = state->getVTherm(Z) * random()->maxwell();
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
    if (state->temperature() > 0. && scatinfo->species < static_cast<int>(2 * numIons))
        lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);

    // Rayleigh scattering in electron rest frame; no support for polarization
    if (scatinfo->species < static_cast<int>(numIons))
    {
        const auto& ion = ions[scatinfo->species];
        _ray->peeloffScattering(I, lambda, ion.Z, pp->direction(), bfkobs);
    }

    // Compton scattering in electron rest frame; with support for polarization if enabled
    else if (scatinfo->species < static_cast<int>(2 * numIons))
    {
        const auto& ion = ions[scatinfo->species - numIons];
        _com->peeloffScattering(I, Q, U, V, lambda, ion.Z, pp->direction(), bfkobs, bfky, pp);
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
    if (state->temperature() > 0. && scatinfo->species < static_cast<int>(2 * numIons))
        lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);

    // room for the outgoing direction
    Direction bfknew;

    // Rayleigh scattering, no support for polarization: determine the new propagation direction
    if (scatinfo->species < static_cast<int>(numIons))
    {
        const auto& ion = ions[scatinfo->species];
        bfknew = _ray->performScattering(lambda, ion.Z, pp->direction());
    }

    // Compton scattering, with support for polarization if enabled:
    // determine the new propagation direction and wavelength, and if polarized, update the stokes vector
    else if (scatinfo->species < static_cast<int>(2 * numIons))
    {
        const auto& ion = ions[scatinfo->species - numIons];
        bfknew = _com->performScattering(lambda, ion.Z, pp->direction(), pp);
    }

    // if we have dispersion, Doppler-shift the outgoing wavelength from the electron rest frame
    if (state->temperature() > 0.) lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfknew, scatinfo->velocity);

    // execute the scattering event in the photon packet
    pp->scatter(bfknew, state->bulkVelocity(), lambda);
}

////////////////////////////////////////////////////////////////////

DisjointWavelengthGrid* XRayIonicGasMix::emissionWavelengthGrid() const
{
    return _emissionGrid;
}

////////////////////////////////////////////////////////////////////

Array XRayIonicGasMix::emissionSpectrum(const MaterialState* state, const Array& /*Jv*/) const
{
    // MAYBE BUILD THE 0 INTO THE TABLE BEFOREHAND?
    Array emis(_numLambda + 2);
    emis[0] = 0.;
    for (int l = 1; l < _numLambda + 1; l++)
    {
        emis[l] = state->getEmissivity(l);
    }
    emis[_numLambda + 1] = 0.;

    return emis * state->volume();
}

// ////////////////////////////////////////////////////////////////////

// Array XRayIonicGasMix::lineEmissionCenters() const
// {
//     return _lambdaL;
// }

// ////////////////////////////////////////////////////////////////////

// Array XRayIonicGasMix::lineEmissionMasses() const
// {
//     return _massL;
// }

// ////////////////////////////////////////////////////////////////////

// Array XRayIonicGasMix::lineEmissionSpectrum(const MaterialState* state, const Array& /*Jv*/) const
// {
//     return _lumL;
// }

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::indicativeTemperature(const MaterialState* state, const Array& /*Jv*/) const
{
    return state->temperature();
}

////////////////////////////////////////////////////////////////////
