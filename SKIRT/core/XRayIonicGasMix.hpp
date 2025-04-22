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

        // can also use file here
        PROPERTY_DOUBLE_LIST(abundances, "the abundances of the ions in the same order as the ions property")

        PROPERTY_DOUBLE(temperature, "the temperature of the gas in K")

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
        IonParam(short Z, short N, double mass, string name) : Z(Z), N(N), mass(mass), name(name) {}

        short Z;      // atomic number
        short N;      // number of electrons
        double mass;  // mass of the ion (amu)
        double vth;   // thermal velocity of the ion (m/s)
        string name;  // name of the ion (eg. Fe1 for neutral iron)
    };

    struct FluorescenceParam
    {
        FluorescenceParam(const Array& a) : Z(a[0]), N(a[1]), n(a[2]), l(a[3]), omega(a[4]), E(a[5]), W(a[6]) {}

        int papIndex{-1};  // index of the photo-absorption transition
        short Z;           // atomic number
        short N;           // number of electrons
        short n;           // principal quantum number of the shell with the hole
        short l;           // orbital quantum number of the subshell with the hole
        double omega;      // fluorescence yield (1)
        double E;          // (central) energy of the emitted photon (eV)
        double W;          // FWHM of the Lorentz shape for the emitted photon (eV), or zero
    };

    struct PhotoAbsorbParam
    {
        PhotoAbsorbParam(const Array& a)
            : Z(a[0]), N(a[1]), n(a[2]), l(a[3]), Eth(a[4]), E0(a[5]), sigma0(a[6]), ya(a[7]), P(a[8]), yw(a[9])
        {}

        int ionIndex{-1};   // index of the ion
        short Z;            // atomic number
        short N;            // number of electrons
        short n;            // principal quantum number of the shell
        short l;            // orbital quantum number of the subshell
        double Eth;         // subshell ionization threshold energy (eV)
        double Emax = 5e5;  // maximum energy for validity of the formula (eV)
        double E0;          // fit parameter (eV)
        double sigma0;      // fit parameter (Mb = 10^-22 m^2)
        double ya;          // fit parameter (1)
        double P;           // fit parameter (1)
        double yw;          // fit parameter (1)
        double y0 = 0.;     // fit parameter (1)
        double y1 = 0.;     // fit parameter (1)

        double Es;
        double sigmamax;

        // return photo-absorption cross section in m2 for given energy in eV and cross section parameters,
        // without taking into account thermal dispersion
        double photoAbsorbSection(double E) const;

        // return photo-absorption cross section in m2 for given energy in eV and cross section parameters,
        // approximating thermal dispersion by replacing the steep threshold transition by a sigmoid
        // error function with given parameters (dispersion and maximum value)
        double photoAbsorbThermalSection(double E) const;
    };

    struct BoundBoundParam
    {
        BoundBoundParam(const Array& a) : Z(a[0]), N(a[1]), E(a[2]), A(a[3]) {}

        int ionIndex{-1};  // index of the ion
        short Z;           // atomic number
        short N;           // number of electrons
        double E;          // energy of the transition (eV)
        double A;          // Einstein A coefficient (s^-1)
        // add Chianti index?
    };

private:
    int _numIons;  // total number of ions
    vector<IonParam> _ionParams;

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

    // thermal velocities and normalized cumulative probability distributions for the scattering channnels:
    //   - Rayleigh scattering by bound electrons for each atom
    //   - Compton scattering by bound electrons for each atom
    //   - fluorescence transitions
    vector<double> _vthermscav;   // indexed on ion + fluo
    ArrayTable<2> _cumprobscavv;  // indexed on wav, 2*ion + fluo

    // bound-electron scattering helpers depending on the configured implementation
    ScatteringHelper* _ray{nullptr};  // Rayleigh scattering helper
    ScatteringHelper* _com{nullptr};  // Compton scattering helper
};

#endif
