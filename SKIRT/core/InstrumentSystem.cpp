/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "InstrumentSystem.hpp"

////////////////////////////////////////////////////////////////////

void InstrumentSystem::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    Instrument* preceding = nullptr;
    for (Instrument* instrument : _instruments)
    {
        if (preceding) instrument->determineSameObserverAsPreceding(preceding);
        preceding = instrument;
    }
}

////////////////////////////////////////////////////////////////////

void InstrumentSystem::flush()
{
    for (Instrument* instrument : _instruments) instrument->flush();
}

////////////////////////////////////////////////////////////////////

void InstrumentSystem::write()
{
    for (Instrument* instrument : _instruments) instrument->write();

    for (RayInstrument* instrument : _rayInstruments) instrument->write();
}

////////////////////////////////////////////////////////////////////

void InstrumentSystem::rayTrace()
{
    for (RayInstrument* instrument: _rayInstruments) instrument->rayTrace();
}