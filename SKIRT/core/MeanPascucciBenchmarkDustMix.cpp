/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MeanPascucciBenchmarkDustMix.hpp"

////////////////////////////////////////////////////////////////////

MaterialMix::ScatteringMode MeanPascucciBenchmarkDustMix::scatteringMode() const
{
    return ScatteringMode::HenyeyGreenstein;
}

//////////////////////////////////////////////////////////////////////

string MeanPascucciBenchmarkDustMix::resourceNameForOpticalProps() const
{
    return "MeanPascucciBenchmarkOpticalProps";
}

//////////////////////////////////////////////////////////////////////