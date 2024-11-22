/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RADIATIONFIELDOPTIONS_HPP
#define RADIATIONFIELDOPTIONS_HPP

#include "DisjointWavelengthGrid.hpp"
#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** The RadiationFieldOptions class simply offers a number of configuration options that are
    related to the radiation field. A simulation always stores the radiation field when it has a
    secondary emission phase or when it has a dynamic medium state (or both). If neither is the
    case, and forced scattering is enabled (see PhotonPacketOptions), the user can still request to
    store the radiation field so that it can be probed for output. */
class RadiationFieldOptions : public SimulationItem
{
    ENUM_DEF(InterpolationMode, NearestNeighbour)
        ENUM_VAL(InterpolationMode, NearestNeighbour,
                 "simplest form of interpolation: take the value of the nearest neighbour")
    ENUM_END()

    ITEM_CONCRETE(RadiationFieldOptions, SimulationItem, "a set of options related to the radiation field")

        PROPERTY_BOOL(storeRadiationField, "store the radiation field so that it can be probed for output")
        ATTRIBUTE_DEFAULT_VALUE(storeRadiationField, "Emission|IteratePrimary:true;false")
        ATTRIBUTE_RELEVANT_IF(storeRadiationField, "!Emission&!IteratePrimary&ForceScattering")
        ATTRIBUTE_DISPLAYED_IF(storeRadiationField, "Level3")
        ATTRIBUTE_INSERT(storeRadiationField,
                         "!Emission&!IteratePrimary&ForceScattering&storeRadiationField:RadiationField")

        PROPERTY_BOOL(storeDirection, "store the direction dependent radiation field to use in reverse ray-tracing")
        ATTRIBUTE_DEFAULT_VALUE(storeDirection, "false")
        ATTRIBUTE_RELEVANT_IF(storeDirection, "storeRadiationField")
        ATTRIBUTE_DISPLAYED_IF(storeDirection, "Level3")

        PROPERTY_ENUM(interpolationMode, InterpolationMode,
                      "what interpolation scheme to use to determine the radiation field in a given direction")
        ATTRIBUTE_DEFAULT_VALUE(interpolationMode, "NearestNeighbour")
        ATTRIBUTE_RELEVANT_IF(interpolationMode, "storeDirection")

        PROPERTY_INT(HEALPIXorder, "HEALPix order: bin amount scales ~ order^2")
        ATTRIBUTE_MIN_VALUE(HEALPIXorder, "0")
        ATTRIBUTE_MAX_VALUE(HEALPIXorder, "15")
        ATTRIBUTE_DEFAULT_VALUE(HEALPIXorder, "6")
        ATTRIBUTE_RELEVANT_IF(HEALPIXorder, "storeDirection")

        PROPERTY_ITEM(radiationFieldWLG, DisjointWavelengthGrid, "the wavelength grid for storing the radiation field")
        ATTRIBUTE_DEFAULT_VALUE(radiationFieldWLG, "LogWavelengthGrid")
        ATTRIBUTE_RELEVANT_IF(radiationFieldWLG, "RadiationField&Panchromatic")

        PROPERTY_BOOL(useReverseRayTracing, "use RRT")
        ATTRIBUTE_DEFAULT_VALUE(useReverseRayTracing, "false")
        ATTRIBUTE_RELEVANT_IF(useReverseRayTracing, "storeDirection")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
