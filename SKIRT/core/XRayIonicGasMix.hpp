/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef XRAYIONICGASMIX_HPP
#define XRAYIONICGASMIX_HPP

#include "ArrayTable.hpp"
#include "EmittingGasMix.hpp"
#include "MaterialMix.hpp"
#include "PhotonPacket.hpp"
#include "SnapshotParameter.hpp"
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

        PROPERTY_STRING(ionNames, "the names of the ions for each element seperated by , (e.g. H1,He2,Fe1,Fe14,...)")

        PROPERTY_DOUBLE_LIST(abundances, "the abundances of the ions in the same order as the ions property")

        PROPERTY_STRING(opticalPropertiesFile, "the name of the file with the optical properties")
        ATTRIBUTE_REQUIRED_IF(opticalPropertiesFile, "True")

        PROPERTY_DOUBLE(temperature, "the temperature of the gas in K")
        ATTRIBUTE_QUANTITY(temperature, "temperature")
        ATTRIBUTE_MIN_VALUE(temperature, "[0")
        ATTRIBUTE_MAX_VALUE(temperature, "1e9]")
        ATTRIBUTE_DEFAULT_VALUE(temperature, "1e4")
        ATTRIBUTE_DISPLAYED_IF(temperature, "Level2")

        PROPERTY_ENUM(scatterBoundElectrons, BoundElectrons, "implementation of scattering by bound electrons")
        ATTRIBUTE_DEFAULT_VALUE(scatterBoundElectrons, "Good")
        ATTRIBUTE_DISPLAYED_IF(scatterBoundElectrons, "Level3")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    explicit XRayIonicGasMix(SimulationItem* parent, string ions, BoundElectrons boundElectrons,
                             vector<double> abundances, double temperature, bool setup);

    void setupSelfBefore() override;

    ~XRayIonicGasMix();

    //======== Private support functions =======

private:
    double interpolateSigma(double lambda, const Array& sigma) const;


    // return thermal velocity for given gas temperature (in K) and particle mass (in amu)
    double vtherm(double amu) const;

    //============= Capabilities =============

    MaterialType materialType() const override;

    bool hasPolarizedScattering() const override;

    bool hasExtraSpecificState() const override;

    bool hasScatteringDispersion() const override;

    bool hasContinuumEmission() const override;

    bool hasLineEmission() const override;

    //============= Medium state setup =============

    vector<StateVariable> specificStateVariableInfo() const override;

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

    void setScatteringInfoIfNeeded(PhotonPacket::ScatteringInfo* scatinfo, double lambda) const;

    void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs, Direction bfky,
                           const MaterialState* state, const PhotonPacket* pp) const override;

    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======== Secondary continuum emission =======

public:
    DisjointWavelengthGrid* emissionWavelengthGrid() const override;

    Array emissionSpectrum(const MaterialState* state, const Array& Jv) const override;

    //======== Secondary line emission =======

public:
    Array lineEmissionCenters() const override;

    Array lineEmissionMasses() const override;

    Array lineEmissionSpectrum(const MaterialState* state, const Array& Jv) const override;

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

    struct Ion
    {
        short Z;  // atomic number
        short N;  // number of electrons

        double mass;   // mass of the ion (amu)
        double vth;    // thermal velocity of the ion (m/s)
        double abund;  // abundance of the ion relative to the total number density
    };

private:
    // all the parameters we want to store, even after the setup
    int _numIons;       // total number of ions (i.e. abundances.size())
    vector<Ion> _ions;  // indexed on ions

    // continuum opacity //
    int _numC;
    Array _lambdaC;
    Array _sigmaabsC;  // indexed on continuum
    Array _sigmascaCO;  // indexed on continuum
    // thermal velocities and normalized cumulative probability distributions for the scattering channnels:
    //   - Rayleigh scattering by bound electrons for each atom
    //   - Compton scattering by bound electrons for each atom
    ArrayTable<2> _cumprobscaCO;  // indexed on continuum, 2*ion

    // continuum emissivity
    DisjointWavelengthGrid* _wavgridCE{nullptr};  // wavelength grid for emission
    Array _emissivityC;                          // indexed on continuum

    // lines //
    Array _lambdaL;
    Array _massL;
    Array _lumL;

    // bound-electron scattering helpers depending on the configured implementation
    ScatteringHelper* _ray{nullptr};  // Rayleigh scattering helper
    ScatteringHelper* _com{nullptr};  // Compton scattering helper
};

#endif
