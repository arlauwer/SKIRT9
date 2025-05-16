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

        PROPERTY_BOOL(levelPopulations, "True if the level populations are used")
        ATTRIBUTE_DEFAULT_VALUE(levelPopulations, "false")

        PROPERTY_STRING(levelPopulationsFilename,
                        "the file containing the level populations of the bound-bound levels in the "
                        "same order as bound-bound levels")
        ATTRIBUTE_RELEVANT_IF(levelPopulationsFilename, "levelPopulations")
        ATTRIBUTE_REQUIRED_IF(levelPopulationsFilename, "levelPopulations")

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
    /** This function returns the index in the private wavelength grid corresponding to the
    specified wavelength. The parameters for converting a wavelength to the appropriate index
    are stored in data members during setup. */
    int indexForLambda(double lambda) const;

    // return thermal velocity for given gas temperature (in K) and particle mass (in amu)
    double vtherm(double amu) const;

    //============= Capabilities =============

    MaterialType materialType() const override;

    bool hasPolarizedScattering() const override;

    bool hasExtraSpecificState() const override;

    bool hasScatteringDispersion() const override;

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

    //======== Secondary emission =======

    /** This function returns a list with the line centers of the supported transitions. */
    Array lineEmissionCenters() const override;

    /** This function returns a list with the masses of the particles emitting each of the lines,
        i.e. the mass of a molecule or atom of the species under consideration in each case. */
    Array lineEmissionMasses() const override;

    /** This function returns a list with the line luminosities for the supported transitions in
        the spatial cell and medium component represented by the specified material state and the
        receiving material mix when it would be embedded in the specified radiation field. */
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

    struct Fluorescence
    {
        double E;  // (central) energy of the emitted photon (eV)
        double W;  // FWHM of the Lorentz shape for the emitted photon (eV), or zero

        short ionIndex;  // index of the ion
    };

    struct BoundTransition
    {
        short upperIndex;  // index of the upper energy level
        short lowerIndex;  // index of the lower energy level
        double lam;        // wavelength of the transition (m)
        double A;          // Einstein A coefficient (s^-1)
        double gul;        // gu/gl

        short ionIndex;  // index of the ion
    };

    struct BoundLevel
    {
        double pop;  // level population relative to parent ion
    };

private:
    // all the parameters we want to store, even after the setup
    int _numIons;       // total number of ions (i.e. abundances.size())
    vector<Ion> _ions;  // indexed on ions

    int _numFluo;                        // total number of fluorescence parameters
    vector<Fluorescence> _fluorescence;  // indexed on photo-absorption parameters

    int _numBLevels;                            // total number of bound-bound levels
    vector<BoundLevel> _boundLevels;            // indexed on bound-bound levels
    int _numBTransitions;                       // total number of bound-bound transitions
    vector<BoundTransition> _boundTransitions;  // indexed on bound-bound transitions

    // all data members are precalculated in setupSelfAfter()

    // wavelength grid (shifted to the left of the actually sampled points to approximate rounding)
    Array _lambdav;  // indexed on wav

    // total extinction and scattering cross sections
    Array _sigmaextv;  // indexed on wav
    Array _sigmascav;  // indexed on wav

    // thermal velocities and normalized cumulative probability distributions for the scattering channnels:
    //   - Rayleigh scattering by bound electrons for each atom
    //   - Compton scattering by bound electrons for each atom
    //   - fluorescence transitions
    ArrayTable<2> _cumprobscavv;  // indexed on wav, 2*ion + fluo

    // bound-electron scattering helpers depending on the configured implementation
    ScatteringHelper* _ray{nullptr};  // Rayleigh scattering helper
    ScatteringHelper* _com{nullptr};  // Compton scattering helper

    TextInFile* _levelPopulationsFile{nullptr};  // file with the level populations of the bound-bound levels
};

#endif
