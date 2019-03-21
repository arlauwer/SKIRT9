/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AbstractWavelengthGridProbe.hpp"

////////////////////////////////////////////////////////////////////

Range AbstractWavelengthGridProbe::wavelengthRange() const
{
    if (wavelengthGrid())
    {
        wavelengthGrid()->setup();
        return Range(wavelengthGrid()->wavelengthRange());
    }
    return Range();
}

////////////////////////////////////////////////////////////////////
