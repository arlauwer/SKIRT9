/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FrameRayInstrument.hpp"
#include "FITSInOut.hpp"
#include "FluxRecorder.hpp"
#include "Log.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PathSegmentGenerator.hpp"
#include "PhotonPacket.hpp"
#include "ProcessManager.hpp"
#include "Table.hpp"

////////////////////////////////////////////////////////////////////

void FrameRayInstrument::setupSelfBefore()
{
    DistantRayInstrument::setupSelfBefore();

    // configure flux recorder
    instrumentRayRecorder()->includeSurfaceBrightnessForDistant(numPixelsX(), numPixelsY(), fieldOfViewX() / numPixelsX(), fieldOfViewY() / numPixelsY(), centerX(), centerY());

    // precalculate information needed by pixelOnDetector() function
    _costheta = cos(inclination());
    _sintheta = sin(inclination());
    _cosphi = cos(azimuth());
    _sinphi = sin(azimuth());
    _cosomega = cos(roll());
    _sinomega = sin(roll());
    _Nxp = numPixelsX();
    _Nyp = numPixelsY();
    _xpmin = centerX() - 0.5 * fieldOfViewX();
    _xpsiz = fieldOfViewX() / numPixelsX();
    _ypmin = centerY() - 0.5 * fieldOfViewY();
    _ypsiz = fieldOfViewY() / numPixelsY();
}

////////////////////////////////////////////////////////////////////

void FrameRayInstrument::rayTrace()
{
    instrumentRayRecorder()->rayTrace();
}

////////////////////////////////////////////////////////////////////