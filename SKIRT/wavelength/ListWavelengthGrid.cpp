/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ListWavelengthGrid.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

ListWavelengthGrid::ListWavelengthGrid(SimulationItem* parent, const vector<double>& wavelengths,
                                       double relativeHalfWidth, bool log, bool setup)
{
    _wavelengths = wavelengths;
    _relativeHalfWidth = relativeHalfWidth;
    _log = log;

    if (setup)
    {
        parent->addChild(this);
        this->setup();
    }
}

////////////////////////////////////////////////////////////////////

void ListWavelengthGrid::setupSelfBefore()
{
    DisjointWavelengthGrid::setupSelfBefore();

    // set the wavelength grid from the list of property values
    if (_relativeHalfWidth)
        setWavelengthBins(NR::array(_wavelengths), _relativeHalfWidth);
    else
        setWavelengthRange(NR::array(_wavelengths), _log);
}

//////////////////////////////////////////////////////////////////////
