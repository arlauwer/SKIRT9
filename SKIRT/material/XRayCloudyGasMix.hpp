/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef XRAYCLOUDYGASMIX_HPP
#define XRAYCLOUDYGASMIX_HPP

#include "CloudyWrapper.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "ElectronScatteringHelper.hpp"
#include "EmittingGasMix.hpp"
#include "ItemInfo.hpp"
#include "MaterialMix.hpp"
#include "PhotonPacket.hpp"

////////////////////////////////////////////////////////////////////

/** The XRayCloudyGasMix class describes the material properties of a partially ionized gas whose
    temperature, ionization state, absorption opacity, continuum emissivity and line emissivities
    are calculated self-consistently by the external photoionization code Cloudy. Each time the
    medium state is updated, a one-zone Cloudy model is evaluated for every spatial cell, given the
    hydrogen density and metallicity of the cell and the radiation field in the cell as calculated
    by SKIRT. Electrons scatter photons through (optionally polarized) Compton scattering.

    <b>Cloudy runs and the lookup table</b>

    The Cloudy runs are performed in subdirectories of the configured table directory. The
    \em tableMode property determines how the results are obtained. In \em Direct mode, Cloudy
    is invoked for every cell and nothing is stored. In \em Record mode, Cloudy is invoked for
    every cell as well, and the results of all runs are stored in a persistent lookup table in the
    table directory; this provides a reference simulation and, at the same time, a table that can
    be used by later simulations or read by external tools. In \em Lookup mode (the default),
    previous runs stored in the table are reused for cells with similar input parameters, both
    within a simulation and across simulations. For each cell, the \em numNeighbors stored runs
    nearest to the cell's input parameters are located. If the nearest run lies within
    \em matchDistance, it is reused as is; otherwise, if all of these runs lie within
    \em interpolationDistance, their results are interpolated; otherwise, Cloudy is invoked and
    the result is added to the table. Refer to the CloudyWrapper class for the definition of the
    distance, the format of the table, and other details. In Record and Lookup mode, the table is
    saved after each medium state update.

    After each medium state update, the log reports how the cells were handled, in the format
    <tt>Cloudy runs: N</tt> (Direct and Record mode) or <tt>Cloudy runs: N, reused: M,
    interpolated: K</tt> (Lookup mode). Cells with zero density are not counted.

    A table directory can only be reused by simulations with the same radiation field wavelength
    grid, minimum radiation field, Cloudy input template, and Cloudy wavelength and line lists.
    For an incompatible configuration, the simulation stops with an error. A table directory
    should not be shared by multiple medium components or simultaneously running simulations.

    Previously saved runs are memory-mapped rather than read into memory, so that large tables
    can be used; runs added during the current simulation are held in memory. When running with
    multiple processes, each process performs the runs for its own cells. After each medium state
    update, the new runs are sent to all processes, and the root process saves them. Within a
    single update, processes may thus perform runs for similar cells that one process would have
    reused or interpolated.

    <b>Resources</b>

    The Cloudy input template, the Cloudy wavelength grid, and the list of emission lines are
    loaded from the resource files <tt>XRayCloudyGasMix_template.in</tt>,
    <tt>XRayCloudyGasMix_wav.dat</tt>, and <tt>XRayCloudyGasMix_lines.dat</tt>. */
class XRayCloudyGasMix : public EmittingGasMix
{
    ENUM_DEF(ElectronScattering, None, Free, FreeWithPolarization)
        ENUM_VAL(ElectronScattering, None, "ignore electron")
        ENUM_VAL(ElectronScattering, Free, "use free-electron Compton scattering for all electrons")
        ENUM_VAL(ElectronScattering, FreeWithPolarization,
                 "use free-electron Compton scattering with support for polarization")
    ENUM_END()

    ENUM_DEF(TableMode, Direct, Record, Lookup)
        ENUM_VAL(TableMode, Direct, "run Cloudy for every cell and store nothing")
        ENUM_VAL(TableMode, Record, "run Cloudy for every cell and store every run in the lookup table")
        ENUM_VAL(TableMode, Lookup, "reuse and interpolate previous Cloudy runs stored in the lookup table")
    ENUM_END()

    ITEM_CONCRETE(XRayCloudyGasMix, EmittingGasMix,
                  "An ionised gas mix using Cloudy to calculate the charge state and temperature self-consistently")
        ATTRIBUTE_TYPE_INSERT(XRayCloudyGasMix, "GasMix,CustomMediumState")

        PROPERTY_ENUM(electronScattering, ElectronScattering, "implementation of scattering by electrons")
        ATTRIBUTE_DEFAULT_VALUE(electronScattering, "Free")
        ATTRIBUTE_DISPLAYED_IF(electronScattering, "Level3")

        PROPERTY_DOUBLE(defaultMetallicity, "default solar metallicity")
        ATTRIBUTE_DEFAULT_VALUE(defaultMetallicity, "1.")
        ATTRIBUTE_MIN_VALUE(defaultMetallicity, "[0")

        PROPERTY_STRING(cloudyExecPath, "path to cloudy executable")
        ATTRIBUTE_DEFAULT_VALUE(cloudyExecPath, "/usr/local/bin/cloudy")

        PROPERTY_STRING(tableDirectory, "path to the cloudy table directory")

        PROPERTY_DOUBLE(radMin, "minimum value of radiation field (W/m2/m)")
        ATTRIBUTE_DEFAULT_VALUE(radMin, "1e-5")
        ATTRIBUTE_MIN_VALUE(radMin, "]0")

        PROPERTY_ENUM(tableMode, TableMode, "how Cloudy runs are obtained and stored")
        ATTRIBUTE_DEFAULT_VALUE(tableMode, "Lookup")
        ATTRIBUTE_DISPLAYED_IF(tableMode, "Level2")

        PROPERTY_INT(numNeighbors, "the number of nearest neighbours used for interpolation")
        ATTRIBUTE_MIN_VALUE(numNeighbors, "1")
        ATTRIBUTE_MAX_VALUE(numNeighbors, "100")
        ATTRIBUTE_DEFAULT_VALUE(numNeighbors, "4")
        ATTRIBUTE_RELEVANT_IF(numNeighbors, "tableModeLookup")
        ATTRIBUTE_DISPLAYED_IF(numNeighbors, "Level3")

        PROPERTY_DOUBLE(matchDistance, "the maximum distance for reusing a previous Cloudy run as is")
        ATTRIBUTE_MIN_VALUE(matchDistance, "[0")
        ATTRIBUTE_DEFAULT_VALUE(matchDistance, "0.01")
        ATTRIBUTE_RELEVANT_IF(matchDistance, "tableModeLookup")
        ATTRIBUTE_DISPLAYED_IF(matchDistance, "Level3")

        PROPERTY_DOUBLE(interpolationDistance, "the maximum distance for interpolating between previous Cloudy runs")
        ATTRIBUTE_MIN_VALUE(interpolationDistance, "[0")
        ATTRIBUTE_DEFAULT_VALUE(interpolationDistance, "0.1")
        ATTRIBUTE_RELEVANT_IF(interpolationDistance, "tableModeLookup")
        ATTRIBUTE_DISPLAYED_IF(interpolationDistance, "Level3")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor releases the scattering helper and the Cloudy wrapper. */
    ~XRayCloudyGasMix();

protected:
    /** This function loads the resource files, sets up the scattering helper, and sets up the
        Cloudy wrapper, loading the lookup table if applicable. */
    void setupSelfBefore() override;

    //======== Private support functions =======

private:
    /** This function returns the index of the Cloudy wavelength bin containing the specified
        wavelength, or -1 if the wavelength lies outside of the Cloudy wavelength grid. */
    int indexForLambda(double lambda) const;

    /** This function returns the cumulative scattering probabilities over the atoms for the
        specified wavelength bin, retrieved from the specified material state. */
    Array getKappaScaCum(const MaterialState* state, int ell) const;

    /** This function stores the cumulative scattering probabilities over the atoms for the
        specified wavelength bin in the specified material state. */
    void setKappaScaCum(MaterialState* state, int ell, const Array& kappaScaCum) const;

    /** This function returns the absorption opacity for the specified wavelength bin, or zero if
        the index is out of range. */
    double kappaAbs(const MaterialState* state, int ell) const;

    /** This function returns the scattering opacity for the specified wavelength bin, or zero if
        the index is out of range or electron scattering is disabled. */
    double kappaSca(const MaterialState* state, int ell) const;

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

    vector<StateVariable> specificStateVariableInfo() const override;

    /** This function initializes the specific state variables of a cell to the values of an empty
        Cloudy result, i.e. all zero. */
    void initializeSpecificState(MaterialState* state, double metallicity, double temperature,
                                 const Array& params) const override;

    //======== Medium state updates =======

    /** This function obtains the Cloudy results for the density, metallicity and radiation field
        of the specified cell, and stores them in the cell's material state. */
    UpdateStatus updateSpecificState(MaterialState* state, const Array& Jv) const override;

    /** This function is called after the medium state of all cells has been updated. It logs
        statistics on the Cloudy runs and, if applicable, shares the new runs between processes
        and saves the lookup table. Because each cell
        update is reported as not converged, the medium state is never considered to be converged,
        so that the number of iterations is determined by the configured maximum. */
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

private:
    void setScatteringInfoIfNeeded(PhotonPacket::ScatteringInfo* scatinfo, double lambda,
                                   const MaterialState* state) const;

public:
    bool peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs, Direction bfky,
                           const MaterialState* state, const PhotonPacket* pp) const override;

    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======== Secondary continuum emission =======

    /** This function returns the Cloudy wavelength grid. */
    DisjointWavelengthGrid* emissionWavelengthGrid() const override;

    /** This function returns the continuum emission spectrum of the specified cell, discretized on
        the extended Cloudy wavelength grid, i.e. with an extra zero value on either side. */
    Array emissionSpectrum(const MaterialState* state, const Array& Jv) const override;

    //======== Secondary line emission =======

    Array lineEmissionCenters() const override;

    Array lineEmissionMasses() const override;

    Array lineEmissionSpectrum(const MaterialState* state, const Array& Jv) const override;

    //======== Temperature =======

    /** This function returns the temperature of the specified cell as calculated by Cloudy. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //======== Private helpers =======

private:
    struct IonParam
    {
        int Z;  // atomic number
        int N;  // number of electrons
    };

    /** This function stores the specified Cloudy results in the specified material state. */
    void storeOutput(MaterialState* state, const Cloudy::Output& output) const;

    //======================== Data Members ========================

private:
    int _numElec{0};  // number of species for electron scattering (one per atom, or none)

    vector<IonParam> _ionParamv;
    Array _lineMassv;
    Array _lineCenterv;

    // MediumState data
    int _indexAbundances{0};       // numIons
    int _indexThermalVelocity{0};  // numAtoms
    int _indexKappaAbs{0};         // numLambda
    int _indexKappaSca{0};         // numLambda (only if _numElec > 0)
    int _indexKappaScaCum{0};      // numLambda x (_numElec-1) (only if _numElec > 0)
    int _indexEmissivity{0};       // numLambda
    int _indexLineEmissivity{0};   // numLines

    CloudyConfig _cloudyConfig;
    CloudyWrapper* _cloudyWrapper{nullptr};

    DisjointWavelengthGrid* _opticalWavelengthGrid{nullptr};

    // compton-electron scattering helpers depending on the configured implementation
    ElectronScatteringHelper::Helper* _com{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
