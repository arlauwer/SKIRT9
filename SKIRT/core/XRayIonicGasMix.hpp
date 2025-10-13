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

////////////////////////////////////////////////////////////////////

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

    ITEM_CONCRETE(XRayIonicGasMix, MaterialMix, "Ionised gas mix")
        ATTRIBUTE_TYPE_INSERT(XRayIonicGasMix, "GasMix,CustomMediumState")

        PROPERTY_STRING(ions, "the names of the ions for each element seperated by , (e.g. H1,He2,Fe1,Fe14,...)")

        PROPERTY_DOUBLE_LIST(abundances, "the abundances of the ions in the same order as the ions property")

        PROPERTY_DOUBLE(temperature, "the temperature of the gas in K")
        ATTRIBUTE_QUANTITY(temperature, "temperature")
        ATTRIBUTE_MIN_VALUE(temperature, "[0")
        ATTRIBUTE_MAX_VALUE(temperature, "1e9]")
        ATTRIBUTE_DEFAULT_VALUE(temperature, "1e4")
        ATTRIBUTE_DISPLAYED_IF(temperature, "Level2")

        PROPERTY_ENUM(scatterBoundElectrons, BoundElectrons, "implementation of scattering by bound electrons")
        ATTRIBUTE_DEFAULT_VALUE(scatterBoundElectrons, "Good")
        ATTRIBUTE_DISPLAYED_IF(scatterBoundElectrons, "Level3")

        PROPERTY_BOOL(resonantScattering, "enable Lyman resonant scattering for all hydrogen-like ions")
        ATTRIBUTE_DEFAULT_VALUE(resonantScattering, "false")
        ATTRIBUTE_DISPLAYED_IF(resonantScattering, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    explicit XRayIonicGasMix(SimulationItem* parent, string ions, vector<double> abundances, double temperature,
                             BoundElectrons boundElectrons, bool setup);

    void setupSelfBefore() override;

    ~XRayIonicGasMix();

    //======== Private support functions =======

private:
    /** This function returns the index in the private wavelength grid corresponding to the
    specified wavelength. The parameters for converting a wavelength to the appropriate index
    are stored in data members during setup. */
    int indexForLambda(double lambda) const;

    //============= Capabilities =============

public:
    MaterialType materialType() const override;

    bool hasPolarizedScattering() const override;

    bool hasResonantScattering() const override;

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

    void setScatteringInfoIfNeeded(PhotonPacket* pp, const MaterialState* state, double lambda) const;

    void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs, Direction bfky,
                           const MaterialState* state, const PhotonPacket* pp) const override;

    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

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

    struct IonParam
    {
        IonParam(short Z, short N) : Z(Z), N(N) {}

        short Z;  // atomic number
        short N;  // number of electrons
    };

private:
    int _numIons;  // total number of ions
    int _numFluo;  // total number of fluorescence transitions
    int _numRes;   // total number of Lyman resonant scattering transitions

    vector<IonParam> _ionParams;  // all the used ion parameters

    // all data members are precalculated in setupSelfAfter()

    // wavelength grid (shifted to the left of the actually sampled points to approximate rounding)
    Array _lambdav;  // indexed on wav

    // total extinction and scattering cross sections
    Array _sigmaextv;  // indexed on wav
    Array _sigmascav;  // indexed on wav

    // emission parameters for each of the fluorescence transitions:
    // if wavelength is nonzero, all photons are emitted at this wavelength;
    // if wavelength is zero, sample wavelength from Lorentz shape defined by central energy and HWHM = FWHM / 2
    vector<double> _lambdafluov;   // indexed on fluo
    vector<double> _centralfluov;  // indexed on fluo
    vector<double> _widthfluov;    // indexed on fluo

    // resonant scattering parameters for each of the Lyman transitions
    vector<double> _aresv;       // indexed on res
    vector<double> _centerresv;  // indexed on res
    vector<int> _Zresv;          // indexed on res

    // thermal velocities and normalized cumulative probability distributions for the scattering channnels:
    //   - Rayleigh scattering by bound electrons for each atom
    //   - Compton scattering by bound electrons for each atom
    //   - fluorescence transitions
    ArrayTable<2> _cumprobscavv;  // indexed on wav, 2*ion + fluo + res

    vector<double> _vtherm;  // indexed on atom

    // bound-electron scattering helpers depending on the configured implementation
    ScatteringHelper* _ray{nullptr};  // Rayleigh scattering helper
    ScatteringHelper* _com{nullptr};  // Compton scattering helper
};

#endif
