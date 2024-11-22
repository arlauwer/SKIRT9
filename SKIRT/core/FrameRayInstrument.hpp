/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FRAMERAYINSTRUMENT_HPP
#define FRAMERAYINSTRUMENT_HPP

#include "DistantRayInstrument.hpp"

////////////////////////////////////////////////////////////////////

/** A FrameRayInstrument object represents a distant instrument that records the surface brightness in
    every pixel of a given frame for each wavelength, and outputs an IFU data cube in a FITS file.

    The instrument allows configuring the field of view and number of pixels in both directions of
    the observer plane. Photon packets arriving from a point that parallel projects outside of the
    field of view are ignored. */
class FrameRayInstrument : public DistantRayInstrument
{
    ITEM_CONCRETE(FrameRayInstrument, DistantRayInstrument,
                  "a distant ray tracing instrument that outputs the surface brightness in every pixel as a data cube")

        PROPERTY_DOUBLE(fieldOfViewX, "the total field of view in the horizontal direction")
        ATTRIBUTE_QUANTITY(fieldOfViewX, "length")
        ATTRIBUTE_MIN_VALUE(fieldOfViewX, "]0")

        PROPERTY_INT(numPixelsX, "the number of pixels in the horizontal direction")
        ATTRIBUTE_MIN_VALUE(numPixelsX, "1")
        ATTRIBUTE_MAX_VALUE(numPixelsX, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsX, "250")

        PROPERTY_DOUBLE(centerX, "the center of the frame in the horizontal direction")
        ATTRIBUTE_QUANTITY(centerX, "length")
        ATTRIBUTE_DEFAULT_VALUE(centerX, "0")
        ATTRIBUTE_DISPLAYED_IF(centerX, "Level2")

        PROPERTY_DOUBLE(fieldOfViewY, "the total field of view in the vertical direction")
        ATTRIBUTE_QUANTITY(fieldOfViewY, "length")
        ATTRIBUTE_MIN_VALUE(fieldOfViewY, "]0")

        PROPERTY_INT(numPixelsY, "the number of pixels in the vertical direction")
        ATTRIBUTE_MIN_VALUE(numPixelsY, "1")
        ATTRIBUTE_MAX_VALUE(numPixelsY, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsY, "250")

        PROPERTY_DOUBLE(centerY, "the center of the frame in the vertical direction")
        ATTRIBUTE_QUANTITY(centerY, "length")
        ATTRIBUTE_DEFAULT_VALUE(centerY, "0")
        ATTRIBUTE_DISPLAYED_IF(centerY, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function configures the FluxRecorder instance associated with this instrument. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

private:

    void rayTrace() override;

    //======================== Data Members ========================

private:
    // data members derived from the discoverable properties during setup, used in pixelOnDetector()
    double _costheta{0};
    double _sintheta{0};
    double _cosphi{0};
    double _sinphi{0};
    double _cosomega{0};
    double _sinomega{0};
    int _Nxp{0};
    int _Nyp{0};
    double _xpmin{0};
    double _xpsiz{0};
    double _ypmin{0};
    double _ypsiz{0};
};

////////////////////////////////////////////////////////////////////

#endif
