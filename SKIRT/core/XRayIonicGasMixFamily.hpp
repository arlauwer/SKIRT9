/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef XRAYIONICGASMIXFAMILY_HPP
#define XRAYIONICGASMIXFAMILY_HPP

#include "XRayIonicGasMix.hpp"
#include "MaterialMixFamily.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the XRayIonicGasMixFamily class represents a family of dust mixes that is
    specified as part of the configuration. Specifically, a property of this class holds a
    user-configurable list of dust mixes representing the family. The family requires a single
    parameter value to select a family member, corresponding to the zero-based index in the
    configured list of dust mixes. The floating point parameter value is rounded to the nearest
    integer and subsequently clipped to be in range. */
class XRayIonicGasMixFamily : public MaterialMixFamily
{
    ENUM_DEF(BoundElectrons, None, Free, FreeWithPolarization, Good, Exact)
        ENUM_VAL(BoundElectrons, None, "ignore bound electrons")
        ENUM_VAL(BoundElectrons, Free, "use free-electron Compton scattering")
        ENUM_VAL(BoundElectrons, FreeWithPolarization,
                 "use free-electron Compton scattering with support for polarization")
        ENUM_VAL(BoundElectrons, Good, "use smooth Rayleigh scattering and exact bound-Compton scattering")
        ENUM_VAL(BoundElectrons, Exact, "use anomalous Rayleigh scattering and exact bound-Compton scattering")
    ENUM_END()
    
    ITEM_CONCRETE(XRayIonicGasMixFamily, MaterialMixFamily, "a family of ionic mixes for each cell")

        PROPERTY_STRING(ions, "the names of the ions for each element seperated by , (e.g. H1,He2,Fe1,Fe14,...)")

        PROPERTY_ENUM(scatterBoundElectrons, BoundElectrons, "implementation of scattering by bound electrons")
        ATTRIBUTE_DEFAULT_VALUE(scatterBoundElectrons, "Good")
        ATTRIBUTE_DISPLAYED_IF(scatterBoundElectrons, "Level3")

    ITEM_END()

    //============= Construction - Setup - Destruction =============
public:

    void setupSelfBefore() override;

    //====================== Other functions =====================

public:
    /** This function returns SnapshotParameter information for a single parameter corresponding to
        the zero-based index in the configured list of dust mixes. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns (a pointer to) the dust mix corresponding to the zero-based index in
        the configured list of dust mixes specified as the single parameter in the list. The
        floating point parameter value is rounded to the nearest integer and subsequently clipped
        to be in range. If the number of parameters in the specified list is not equal to one, the
        behavior is undefined. The material mix family retains ownership of the returned dust mix,
        and guarantees that it will not be destroyed until the family itself is destroyed. */
    const MaterialMix* mix(const Array& parameters) override;

    //======================== Data Members ========================
private:
    vector<string> _ionNames;

};

////////////////////////////////////////////////////////////////////

#endif
