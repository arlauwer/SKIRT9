/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef XRAYIONICGASMIX_HPP
#define XRAYIONICGASMIX_HPP

#include "ArrayTable.hpp"
#include "MaterialMix.hpp"
#include "PhotonPacket.hpp"
#include "SnapshotParameter.hpp"

////////////////////////////////////////////////////////////////////

struct IonParams;
struct FluorescenceParams;

class XRayIonicGasMix : public MaterialMix
{
    ENUM_DEF(BoundElectrons, None, Free, FreeWithPolarization, Good, Exact)
        ENUM_VAL(BoundElectrons, None, "ignore bound electrons")
        ENUM_VAL(BoundElectrons, Free, "use free-electron Compton scattering")
        ENUM_VAL(BoundElectrons, FreeWithPolarization,
                 "use free-electron Compton scattering with support for polarization")
        ENUM_VAL(BoundElectrons, Good, "use smooth Rayleigh scattering and exact bound-Compton scattering")
        ENUM_VAL(BoundElectrons, Exact, "use anomalous Rayleigh scattering and exact bound-Compton scattering")
    ENUM_END()

    ITEM_CONCRETE(XRayIonicGasMix, MaterialMix,
                  "Ionised gas mix")

        PROPERTY_STRING(filename, "the name of the file with ion abundancies per zone")

        PROPERTY_STRING(ions, "the names of the ions for each element seperated by , (e.g. H1,He2,Fe1,Fe14,...)")

        PROPERTY_DOUBLE_LIST(abundancies, "the abundancies for the ions")
        ATTRIBUTE_MIN_VALUE(abundancies, "[0")
        ATTRIBUTE_MAX_VALUE(abundancies, "1]")
        ATTRIBUTE_REQUIRED_IF(abundancies, "false")
        ATTRIBUTE_DISPLAYED_IF(abundancies, "Level2")

        PROPERTY_ENUM(scatterBoundElectrons, BoundElectrons, "implementation of scattering by bound electrons")
        ATTRIBUTE_DEFAULT_VALUE(scatterBoundElectrons, "Good")
        ATTRIBUTE_DISPLAYED_IF(scatterBoundElectrons, "Level3")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    void setupSelfBefore() override;

    ~XRayIonicGasMix();

    //============= Capabilities =============

    MaterialType materialType() const override;

    bool hasPolarizedScattering() const override;

    bool hasExtraSpecificState() const override;

    bool hasScatteringDispersion() const override;

    //============= Medium state setup =============

    vector<SnapshotParameter> parameterInfo() const override;

    vector<StateVariable> specificStateVariableInfo() const override;

    void initializeSpecificState(MaterialState* state, double metallicity, double temperature,
                                 const Array& params) const override;

    //============= Low-level material properties =============

    double mass() const override;

    double sectionAbs(double lambda) const override;

    double sectionSca(double lambda) const override;

    double sectionExt(double lambda) const override;

    //============= High-level photon life cycle =============

    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    Array sigmaSca(double lambda, const MaterialState* state) const;

    void setScatteringInfoIfNeeded(PhotonPacket::ScatteringInfo* scatinfo, double lambda,
                                   const MaterialState* state) const;

    void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs, Direction bfky,
                           const MaterialState* state, const PhotonPacket* pp) const override;

    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======================== Data Members ========================

public:
    // base class for bound-electron scattering helpers (public because we derive from it in anonymous namespace)
    class ScatteringHelper;

private:
    // all data members are precalculated in setupSelfAfter()

    int _numIons;   // total number of ions
    int _numAtoms;  // all unique Z
    int _numFluos;  // number of fluorescence transitions

    Array _fluorescenceMass;  // mass of ion for each fluorescence transition, used for setScatteringInfoIfNeeded

    vector<IonParams*> _ionParams;
    vector<const FluorescenceParams*> _fluorescenceParams;  // used for setScatteringInfoIfNeeded

    // bound-electron scattering helpers depending on the configured implementation
    ScatteringHelper* _ray{nullptr};  // Rayleigh scattering helper
    ScatteringHelper* _com{nullptr};  // Compton scattering helper
};

#endif
