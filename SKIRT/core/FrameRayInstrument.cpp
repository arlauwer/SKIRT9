/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FrameRayInstrument.hpp"
#include "FluxRecorder.hpp"
#include "PhotonPacket.hpp"
#include "Table.hpp"
#include "Log.hpp"
#include "ParallelFactory.hpp"
#include "Parallel.hpp"

////////////////////////////////////////////////////////////////////

void FrameRayInstrument::setupSelfBefore()
{
    DistantRayInstrument::setupSelfBefore();

    // configure flux recorder
    // instrumentFluxRecorder()->includeSurfaceBrightnessForDistant(numPixelsX(), numPixelsY(), fieldOfViewX() / numPixelsX(), fieldOfViewY() / numPixelsY(), centerX(), centerY());

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

int FrameRayInstrument::pixelOnDetector(const PhotonPacket* pp) const
{
    // get the position
    double x, y, z;
    pp->position().cartesian(x, y, z);

    // transform to detector coordinates using inclination, azimuth, and roll angle
    double xpp = -_sinphi * x + _cosphi * y;
    double ypp = -_cosphi * _costheta * x - _sinphi * _costheta * y + _sintheta * z;
    double xp = _cosomega * xpp - _sinomega * ypp;
    double yp = _sinomega * xpp + _cosomega * ypp;

    // scale and round to pixel index
    int i = static_cast<int>(floor((xp - _xpmin) / _xpsiz));
    int j = static_cast<int>(floor((yp - _ypmin) / _ypsiz));
    if (i < 0 || i >= _Nxp || j < 0 || j >= _Nyp)
        return -1;
    else
        return i + _Nxp * j;
}

////////////////////////////////////////////////////////////////////

void FrameRayInstrument::rayTrace()
{
    // get frame configuration
    int Nxp = numPixelsX();
    int Nyp = numPixelsY();
    double xcent = centerX();
    double ycent = centerY();
    double xpmin = centerX() - 0.5 * fieldOfViewX();
    double xpsiz = fieldOfViewX() / numPixelsX();
    double ypmin = centerY() - 0.5 * fieldOfViewY();
    double ypsiz = fieldOfViewY() / numPixelsY();

    // calculate sines and cosines for the transformation from observer to model coordinates
    double costheta = cos(inclination());
    double sintheta = sin(inclination());
    double cosphi = cos(azimuth());
    double sinphi = sin(azimuth());
    double cosomega = cos(roll());
    double sinomega = sin(roll());

    // determine the observer axes as directions in model coordinates (inverse transform of frame instrument)
    // k_z is the direction from observer to model (i.e. left-handed coordinate frame)
    Direction kz(-cosphi * sintheta, -sinphi * sintheta, -costheta, false);
    Direction ky(-cosphi * costheta * cosomega - sinphi * sinomega, -sinphi * costheta * cosomega + cosphi * sinomega,
                 sintheta * cosomega, false);
    Direction kx(cosphi * costheta * sinomega - sinphi * cosomega, sinphi * costheta * sinomega + cosphi * cosomega,
                 -sintheta * sinomega, false);

    // get a distance that should be well outside of the model (but with a similar order of magnitude)
    double zp = 10. * (fieldOfViewX() + fieldOfViewY());

    // allocate result array with the appropriate size and initialize contents to zero
    Table<2> vvv(Nyp, Nxp);  // reverse index order to get proper data value ordering for FITSInOut::write()

    // calculate the results in parallel
    auto log = find<Log>();
    log->infoSetElapsed(Nyp);
    auto parallel = find<ParallelFactory>()->parallelDistributed();
    parallel->call(Nyp, [&vvv, xpmin, xpsiz, ypmin, ypsiz, kx, ky, kz, zp, costheta, sintheta, cosphi, sinphi,
                         cosomega, sinomega, Nxp, log](size_t firstIndex, size_t numIndices) {
        double val;
        string progress = "calculated projected pixels: ";

        // loop over pixels
        for (size_t j = firstIndex; j != firstIndex + numIndices; ++j)
        {
            for (int i = 0; i != Nxp; ++i)
            {
                // transform pixel indices to observer coordinates
                double xp = xpmin + (i + 0.5) * xpsiz;
                double yp = ypmin + (j + 0.5) * ypsiz;

                // transform observer coordinates to model coordinates (inverse transform of frame instrument)
                double xpp = sinomega * xp - cosomega * yp;
                double ypp = cosomega * xp + sinomega * yp;
                double zpp = zp;
                double x = cosphi * costheta * xpp - sinphi * ypp + cosphi * sintheta * zpp;
                double y = sinphi * costheta * xpp + cosphi * ypp + sinphi * sintheta * zpp;
                double z = -sintheta * xpp + costheta * zpp;

                // get the quantity value aggregated along a path leaving the pixel through the model
                // bridge->valuesAlongPath(Position(x, y, z), kz, values);


            }
            log->infoIfElapsed(progress, 1);
        }
    });
}

////////////////////////////////////////////////////////////////////