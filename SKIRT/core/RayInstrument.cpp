/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "RayInstrument.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "FluxRecorder.hpp"

////////////////////////////////////////////////////////////////////

void RayInstrument::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // select "local" or default wavelength grid
    auto config = find<Configuration>();
    _instrumentWavelengthGrid = config->wavelengthGrid(wavelengthGrid());

    // discover details about the simulation
    bool hasMedium = config->hasMedium();
    bool hasMediumEmission = config->hasSecondaryEmission();
}

////////////////////////////////////////////////////////////////////

std::string RayInstrument::itemName() const
{
    return instrumentName();
}

////////////////////////////////////////////////////////////////////

void RayInstrument::write()
{
    
}

////////////////////////////////////////////////////////////////////
