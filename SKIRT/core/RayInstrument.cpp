/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "RayInstrument.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "FluxRecorder.hpp"
#include "PathSegmentGenerator.hpp"
#include "SpatialGridPath.hpp"

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
    _ms = find<MediumSystem>();
    _grid = _ms->grid();

    _recorder = new RayRecorder(this);
    _recorder->setSimulationInfo(instrumentName(), instrumentWavelengthGrid(), hasMedium, hasMediumEmission);
    _recorder->setUserFlags(_recordComponents, _numScatteringLevels, _recordPolarization, _recordStatistics);
}

RayInstrument::~RayInstrument()
{
    delete _recorder;
}

////////////////////////////////////////////////////////////////////

void RayInstrument::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    // finalize configuration of the flux recorder
    _recorder->finalizeConfiguration();
}

////////////////////////////////////////////////////////////////////


std::string RayInstrument::itemName() const
{
    return instrumentName();
}

////////////////////////////////////////////////////////////////////

void RayInstrument::write()
{
    _recorder->calibrateAndWrite();
}

////////////////////////////////////////////////////////////////////
