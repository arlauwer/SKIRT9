/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LyaNeutralHydrogenGasMix.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "LyUtils.hpp"
#include "MaterialState.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // the combined Lya1 and Lya2 for hydrogen
    constexpr int Z = 1;
    constexpr double lyaA = Constants::EinsteinALya();
    constexpr double center = Constants::lambdaLya();
    constexpr double g = 3.;
    constexpr double m = Constants::Mproton();
}

////////////////////////////////////////////////////////////////////

void LyaNeutralHydrogenGasMix::setupSelfBefore()
{
    MaterialMix::setupSelfBefore();

    _dpf.initialize(random(), includePolarization());
    _vth = LyUtils::vtherm(defaultTemperature(), m);
}

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType LyaNeutralHydrogenGasMix::materialType() const
{
    return MaterialType::Gas;
}

////////////////////////////////////////////////////////////////////

bool LyaNeutralHydrogenGasMix::hasPolarizedScattering() const
{
    return includePolarization();
}

////////////////////////////////////////////////////////////////////

bool LyaNeutralHydrogenGasMix::hasResonantScattering() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool LyaNeutralHydrogenGasMix::hasExtraSpecificState() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool LyaNeutralHydrogenGasMix::hasScatteringDispersion() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> LyaNeutralHydrogenGasMix::specificStateVariableInfo() const
{
    vector<StateVariable> result = {StateVariable::numberDensity(), StateVariable::temperature()};
    result.push_back(StateVariable::custom(0, "thermal velocity", "velocity"));
    result.push_back(StateVariable::custom(0, "voigt parameter", ""));
    return result;
}

////////////////////////////////////////////////////////////////////

void LyaNeutralHydrogenGasMix::initializeSpecificState(MaterialState* state, double /*metallicity*/, double temperature,
                                                       const Array& /*params*/) const
{
    // leave the temperature at zero if the cell does not contain any material for this component
    if (state->numberDensity() > 0.)
    {
        // if no temperature was imported, use default value
        if (temperature < 0) temperature = defaultTemperature();

        // make sure the temperature is at least the local universe CMB temperature
        state->setTemperature(max(Constants::Tcmb(), temperature));

        double vth = LyUtils::vtherm(temperature, m);
        double a = lyaA * center / 4. / M_PI / vth;
        state->setCustom(0, vth);
        state->setCustom(1, a);
    }
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenGasMix::mass() const
{
    return Constants::Mproton();
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenGasMix::sectionAbs(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenGasMix::sectionSca(double lambda) const
{
    double a = lyaA * center / 4. / M_PI / _vth;
    return LyUtils::section(_vth, a, center, g, lambda);
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenGasMix::sectionExt(double lambda) const
{
    double a = lyaA * center / 4. / M_PI / _vth;
    return LyUtils::section(_vth, a, center, g, lambda);
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenGasMix::opacityAbs(double /*lambda*/, const MaterialState* /*state*/,
                                            const PhotonPacket* /*pp*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenGasMix::opacitySca(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double n = state->numberDensity();
    double vth = state->custom(0);
    double a = state->custom(1);
    return n > 0. ? n * LyUtils::section(vth, a, center, g, lambda) : 0.;
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenGasMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double n = state->numberDensity();
    double vth = state->custom(0);
    double a = state->custom(1);
    return n > 0. ? n * LyUtils::section(vth, a, center, g, lambda) : 0.;
}

////////////////////////////////////////////////////////////////////

void LyaNeutralHydrogenGasMix::peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda,
                                                 Direction bfkobs, Direction bfky, const MaterialState* state,
                                                 const PhotonPacket* pp) const
{
    // draw a random atom velocity & phase function, unless a previous peel-off stored this already
    auto scatinfo = const_cast<PhotonPacket*>(pp)->getScatteringInfo();
    if (!scatinfo->valid)
    {
        scatinfo->valid = true;
        double vth = state->custom(0);
        double a = state->custom(1);
        std::tie(scatinfo->velocity, scatinfo->dipole) = LyUtils::sampleAtomVelocity(
            vth, a, center, lambda, state->temperature(), state->numberDensity(), pp->direction(), config(), random());
    }

    // add the contribution to the Stokes vector components depending on scattering type
    if (scatinfo->dipole)
    {
        // contribution of dipole scattering with or without polarization
        _dpf.peeloffScattering(I, Q, U, V, pp->direction(), bfkobs, bfky, pp);
    }
    else
    {
        // isotropic scattering removes polarization, so the contribution is trivially 1
        I += 1.;
    }

    // Doppler-shift the photon packet wavelength into and out of the atom frame
    lambda = LyUtils::shiftWavelength(lambda, scatinfo->velocity, pp->direction(), bfkobs);
}

////////////////////////////////////////////////////////////////////

void LyaNeutralHydrogenGasMix::performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const
{
    // draw a random atom velocity & phase function, unless a peel-off stored this already
    auto scatinfo = pp->getScatteringInfo();
    if (!scatinfo->valid)
    {
        scatinfo->valid = true;
        double vth = state->custom(0);
        double a = state->custom(1);
        std::tie(scatinfo->velocity, scatinfo->dipole) = LyUtils::sampleAtomVelocity(
            vth, a, center, lambda, state->temperature(), state->numberDensity(), pp->direction(), config(), random());
    }

    // draw the outgoing direction from the dipole or the isotropic phase function
    // and, if required, update the polarization state of the photon packet
    Direction bfknew;
    if (scatinfo->dipole)
    {
        bfknew = _dpf.performScattering(pp->direction(), pp);
    }
    else
    {
        bfknew = random()->direction();
        if (includePolarization()) pp->setUnpolarized();
    }

    // Doppler-shift the photon packet wavelength into and out of the atom frame
    lambda = LyUtils::shiftWavelength(lambda, scatinfo->velocity, pp->direction(), bfknew);

    // execute the scattering event in the photon packet
    pp->scatter(bfknew, state->bulkVelocity(), lambda);
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenGasMix::indicativeTemperature(const MaterialState* state, const Array& /*Jv*/) const
{
    return state->temperature();
}

////////////////////////////////////////////////////////////////////
