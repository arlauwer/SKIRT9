/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef XRAYIONICGASMIX_HPP
#define XRAYIONICGASMIX_HPP

#include "CloudyWrapper.hpp"
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

        PROPERTY_STRING(libLocation, "Cloudy library location")
        ATTRIBUTE_DISPLAYED_IF(libLocation, "Level3")

        PROPERTY_ENUM(scatterBoundElectrons, BoundElectrons, "implementation of scattering by bound electrons")
        ATTRIBUTE_DEFAULT_VALUE(scatterBoundElectrons, "Good")
        ATTRIBUTE_DISPLAYED_IF(scatterBoundElectrons, "Level2")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    void setupSelfBefore() override;

    ~XRayIonicGasMix();

    //======== Private support functions =======

private:
    int indexForOpacity(double lambda) const;

    //============= Capabilities =============

public:
    MaterialType materialType() const override;

    bool hasPolarizedScattering() const override;

    bool hasExtraSpecificState() const override;

    DynamicStateType hasDynamicMediumState() const override;

    bool hasScatteringDispersion() const override;

    bool hasContinuumEmission() const override;

    bool hasLineEmission() const override;

    //======== Medium state setup =======

public:
    vector<StateVariable> specificStateVariableInfo() const override;

    void initializeSpecificState(MaterialState* state, double metallicity, double temperature,
                                 const Array& params) const override;

    //======== Medium state updates =======

    UpdateStatus updateSpecificState(MaterialState* state, const Array& Jv) const override;

    bool isSpecificStateConverged(int numCells, int numUpdated, int numNotConverged, MaterialState* currentAggregate,
                                  MaterialState* previousAggregate) const override;

private:
    double updateState(MaterialState* state, double n, double Z, const Array& J);

    //============= Low-level material properties =============

public:
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

    void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs, Direction bfky,
                           const MaterialState* state, const PhotonPacket* pp) const override;

    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======== Secondary continuum emission =======

public:
    DisjointWavelengthGrid* emissionWavelengthGrid() const override;

    Array emissionSpectrum(const MaterialState* state, const Array& Jv) const override;

    //======== Secondary line emission =======

public:
    // Array lineEmissionCenters() const override;

    // Array lineEmissionMasses() const override;

    // Array lineEmissionSpectrum(const MaterialState* state, const Array& Jv) const override;

    //======== Temperature =======

public:
    /** This function returns an indicative temperature of the material mix when it would be
    embedded in a given radiation field. The implementation in this class ignores the radiation
    field and returns the (spatially constant) temperature configured for this material mix. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //======================== Data Members ========================

public:
    // base class for bound-electron scattering helpers (public because we derive from it in anonymous namespace)
    class ScatteringHelper;

private:
    Log* _log{nullptr};

    Range _range;

    // MaterialState indices: size //
    // Abundances:      numIons
    // VTherm:          numAtoms
    // SigmaAbs:        numLambda
    // SigmaSca:        numLambda
    // SigmaScaCum:     numLambda x 2*numIons+1 (+1 for cumulative)
    // Emissivity:      numLambda + 2 (+2 for 'outside' bins)
    int _indexAbundances;       // index of the abundances in the custom state variables
    int _indexThermalVelocity;  // index of the thermal velocity in the custom state variables
    int _indexKappaAbs;         // index of the absorption cross section in the custom state variables
    int _indexKappaSca;         // index of the absorption and scattering cross sections in the custom state variables
    int _indexKappaScaCum;      // index of the cumulative scattering cross section for Rayleigh and Compton scattering
    int _indexEmissivity;       // index of the emissivity in the custom state variables

    // shared variables //
    Array _lambda;

    // int _numLambdaEmis;
    // Array _lambdaEmis;
    DisjointWavelengthGrid* _emissionGrid{nullptr};  // grid for emission

    // bound-electron scattering helpers depending on the configured implementation
    ScatteringHelper* _ray{nullptr};  // Rayleigh scattering helper
    ScatteringHelper* _com{nullptr};  // Compton scattering helper

    // INFO: the density (1/m3) is the hden in Cloudy!
    // Meaning that if you set the density in SKIRT (mix.txt) to 1e4, it will be the hden in Cloudy!
    // the sum of the H abundances: H, H+ (H+ not present in SKIRT yet?), (molecules currently turned off) will be 1e4

    CloudyWrapper _cloudyWrapper;
};

#endif
