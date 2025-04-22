/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SelectDustMixFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> SelectDustMixFamily::parameterInfo() const
{
    return {SnapshotParameter::custom("dustmix index")};
}

////////////////////////////////////////////////////////////////////

const MaterialMix* SelectDustMixFamily::mix(double /*Z*/, double /*T*/, const Array& parameters)
{
    long numMixes = _dustMixes.size();
    long index = max(0L, min(std::lround(parameters[0]), numMixes - 1));
    return _dustMixes[index];
}

////////////////////////////////////////////////////////////////////

const MaterialMix* SelectDustMixFamily::mix()
{
    return _dustMixes[0];
}

////////////////////////////////////////////////////////////////////