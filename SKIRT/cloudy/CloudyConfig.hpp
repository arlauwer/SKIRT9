#ifndef CLOUDY_CONFIG_HPP
#define CLOUDY_CONFIG_HPP

#include "Array.hpp"
#include "Basics.hpp"
#include "DisjointWavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

struct CloudyConfig
{
    struct RadField
    {
        int numBins{0};
        Array edgev;         // rydberg, ascending, size numBins+1
        Array widthv;        // m, size numBins
        double minFlux{0.};  // W/m2/m, floor used to avoid zeros in the SED file
    };

    struct WavGrid
    {
        int numBins{0};
        Array lambdav;  // m, ascending, size numBins
    };

    RadField rad;
    WavGrid wav;

    int numIons{0};   // number of ion slots in the SKIRT abundance vector
    int numLines{0};  // number of lines requested by the template

    string execPath;  // path to the cloudy executable

    void setup(const DisjointWavelengthGrid& radGrid, const DisjointWavelengthGrid& wavGrid, double radMin, int numIons,
               int numLines, string execPath);
};

////////////////////////////////////////////////////////////////////

#endif
