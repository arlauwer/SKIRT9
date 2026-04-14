/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef XRAYIONICGASMIX_HPP
#define XRAYIONICGASMIX_HPP

#include "CloudyWrapper.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "EmittingGasMix.hpp"
#include "ItemInfo.hpp"
#include "Log.hpp"
#include "MaterialMix.hpp"
#include "PhotonPacket.hpp"
#include "TextInFile.hpp"

////////////////////////////////////////////////////////////////////

class XRayIonicGasMix : public EmittingGasMix
{
    ENUM_DEF(BoundElectrons, None, Free, FreeWithPolarization, Good, Exact)
        ENUM_VAL(BoundElectrons, None, "ignore bound electrons")
        ENUM_VAL(BoundElectrons, Free, "use free-electron Compton scattering")
        ENUM_VAL(BoundElectrons, FreeWithPolarization,
                 "use free-electron Compton scattering with support for polarization")
        ENUM_VAL(BoundElectrons, Good, "use smooth Rayleigh scattering and exact bound-Compton scattering")
        ENUM_VAL(BoundElectrons, Exact, "use anomalous Rayleigh scattering and exact bound-Compton scattering")
    ENUM_END()

    ITEM_CONCRETE(XRayIonicGasMix, EmittingGasMix, "Ionised gas mix")
        ATTRIBUTE_TYPE_INSERT(XRayIonicGasMix, "GasMix,CustomMediumState")

        PROPERTY_ENUM(scatterBoundElectrons, BoundElectrons, "implementation of scattering by bound electrons")
        ATTRIBUTE_DEFAULT_VALUE(scatterBoundElectrons, "Good")
        ATTRIBUTE_DISPLAYED_IF(scatterBoundElectrons, "Level2")

        PROPERTY_ITEM(opticalWavelengthGrid, DisjointWavelengthGrid, "optical wavelength grid")

        PROPERTY_STRING(cloudyExecPath, "path to cloudy executable")
        ATTRIBUTE_DEFAULT_VALUE(cloudyExecPath, "/usr/local/bin/cloudy")

        PROPERTY_DOUBLE(radMin, "minimum value of radiation field (W/m2/m)")
        ATTRIBUTE_DEFAULT_VALUE(radMin, "1e-5")
        ATTRIBUTE_MIN_VALUE(radMin, "]0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    void setupSelfBefore() override;

    ~XRayIonicGasMix();

    //============= Capabilities =============

    MaterialType materialType() const override;

    bool hasPolarizedScattering() const override;

    bool hasExtraSpecificState() const override;

    DynamicStateType hasDynamicMediumState() const override;

    bool hasScatteringDispersion() const override;

    bool hasContinuumEmission() const override;

    bool hasLineEmission() const override;

    //======== Medium state setup =======

    vector<StateVariable> specificStateVariableInfo() const override;

    void initializeSpecificState(MaterialState* state, double metallicity, double temperature,
                                 const Array& params) const override;

    //======== Medium state updates =======

    UpdateStatus updateSpecificState(MaterialState* state, const Array& Jv) const override;

    bool isSpecificStateConverged(int numCells, int numUpdated, int numNotConverged, MaterialState* currentAggregate,
                                  MaterialState* previousAggregate) const override;

    //============= Low-level material properties =============

    double mass() const override;

    double sectionAbs(double lambda) const override;

    double sectionSca(double lambda) const override;

    double sectionExt(double lambda) const override;

    //============= High-level photon life cycle =============

    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    void setScatteringInfoIfNeeded(PhotonPacket::ScatteringInfo* scatinfo, double lambda,
                                   const MaterialState* state) const;

    bool peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs, Direction bfky,
                           const MaterialState* state, const PhotonPacket* pp) const override;

    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======== Secondary continuum emission =======

    DisjointWavelengthGrid* emissionWavelengthGrid() const override;

    Array emissionSpectrum(const MaterialState* state, const Array& Jv) const override;

    //======== Secondary line emission =======

    Array lineEmissionCenters() const override;

    Array lineEmissionMasses() const override;

    Array lineEmissionSpectrum(const MaterialState* state, const Array& Jv) const override;

    //======== Temperature =======

    /** This function returns an indicative temperature of the material mix when it would be
    embedded in a given radiation field. The implementation in this class ignores the radiation
    field and returns the (spatially constant) temperature configured for this material mix. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //======== Helper Functions =======

private:
    void setupCloudyConfig();

    void updateSpecificState(MaterialState* state, const Cloudy::Output& output) const;

    //======================== Data Members ========================

public:
    // base class for bound-electron scattering helpers (public because we derive from it in anonymous namespace)
    class ScatteringHelper;

private:
    Log* _log{nullptr};

    // bound-electron scattering helpers depending on the configured implementation
    ScatteringHelper* _ray{nullptr};  // Rayleigh scattering helper
    ScatteringHelper* _com{nullptr};  // Compton scattering helper

    // MediumState data
    int _indexAbundances;       // numIons
    int _indexThermalVelocity;  // numAtoms
    int _indexKappaAbs;         // numLambda
    int _indexKappaSca;         // numLambda
    int _indexKappaScaCum;      // numLambda x 2*numIons+1 (+1 for cumulative)
    int _indexEmissivity;       // numLambda + 2 (+2 for 'outside' bins)
    int _indexLineEmissivity;   // numLines

    CloudyConfig _cloudyConfig;
    CloudyWrapper _cloudyWrapper;

    DisjointWavelengthGrid* _emissionWavelengthGrid;
};

#endif
