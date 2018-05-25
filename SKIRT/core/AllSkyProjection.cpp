/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AllSkyProjection.hpp"

////////////////////////////////////////////////////////////////////

void AllSkyProjection::fromSphereToRectangle(double inclination, double azimuth, double& x, double& y) const
{
    double longitude = -azimuth;
    double latitude = inclination - M_PI_2;
    fromGlobeToRectangle(longitude, latitude, x, y);
}

////////////////////////////////////////////////////////////////////

bool AllSkyProjection::fromRectangleToSphere(double x, double y, double& inclination, double& azimuth) const
    double longitude, latitude;
    bool success = fromRectangleToGlobe(x, y, longitude, latitude);
    inclination = latitude + M_PI_2;
    azimuth = -longitude;
    return success;
}

////////////////////////////////////////////////////////////////////
