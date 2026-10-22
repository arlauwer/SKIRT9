/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CLOUDYCONFIG_HPP
#define CLOUDYCONFIG_HPP

#include "Array.hpp"
#include "Basics.hpp"
#include <functional>

////////////////////////////////////////////////////////////////////

/** CloudyConfig holds the configuration shared by all Cloudy runs of a single client, i.e. the
    discretization of the input radiation field, the wavelength grid on which Cloudy reports
    opacities and emissivities, the number of ions and lines in the output, and the location of
    the Cloudy executable. The client fills in all fields before handing the configuration to a
    CloudyWrapper instance. */
struct CloudyConfig
{
    /** The discretization of the radiation field passed to Cloudy. */
    struct RadField
    {
        int numBins{0};
        Array borderv;       // Ry, size numBins+1, same order as the bins (i.e. descending energy)
        Array widthv;        // m, size numBins
        double minFlux{0.};  // W/m2/m/sr, floor used to avoid zeros in the SED file
    };

    /** The wavelength grid on which Cloudy reports opacities and emissivities, which must match
        the energy grid used internally by Cloudy. */
    struct WavGrid
    {
        int numBins{0};
        Array lambdav;  // m, ascending, size numBins
    };

    RadField rad;
    WavGrid wav;

    int numIons{0};   // number of ion slots in the abundance vector
    int numLines{0};  // number of lines requested by the template

    /** Returns the index in the abundance vector for the ion with the given Cloudy species name,
        or -1 if the species should be ignored. */
    std::function<int(const string&)> ionIndex;

    string execPath;  // path to the Cloudy executable
};

////////////////////////////////////////////////////////////////////

#endif
