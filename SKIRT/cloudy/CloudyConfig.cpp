#include "CloudyConfig.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

void CloudyConfig::setup(const DisjointWavelengthGrid& radGrid, const DisjointWavelengthGrid& wavGrid, double radMin,
                         int numIons, int numLines, string execPath)
{
    constexpr double front = Constants::h() * Constants::c() / Constants::Qelectron();

    rad.numBins = radGrid.numBins();
    rad.edgev = front / radGrid.borderv();
    rad.widthv = radGrid.dlambdav();
    rad.minFlux = radMin;

    wav.numBins = wavGrid.numBins();
    wav.lambdav = wavGrid.lambdav();

    this->numIons = numIons;
    this->numLines = numLines;

    this->execPath = execPath;
}

////////////////////////////////////////////////////////////////////
