/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "XRayCloudyGasMix.hpp"
#include "AtomUtils.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "ListBorderWavelengthGrid.hpp"
#include "Log.hpp"
#include "MaterialState.hpp"
#include "NR.hpp"
#include "ProcessManager.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "TextInFile.hpp"
#include "UpdateStatus.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    constexpr int numAtoms = 30;
    constexpr int numIons = 495;

    // returns the thermal velocity for the given temperature and particle mass
    double vtherm(double T, double mass)
    {
        return sqrt(Constants::k() * T / mass);
    }
}

////////////////////////////////////////////////////////////////////

void XRayCloudyGasMix::setupSelfBefore()
{
    EmittingGasMix::setupSelfBefore();

    // setup all ions
    _ionParamv.resize(numIons);
    for (int Z = 1; Z <= numAtoms; Z++)
    {
        for (int N = 0; N <= Z; N++)
        {
            int i = AtomUtils::ionIndex(Z, N);
            _ionParamv[i].Z = Z;
            _ionParamv[i].N = N;
        }
    }

    // create scattering helpers depending on the user-configured implementation type
    using namespace ElectronScatteringHelper;
    switch (electronScattering())
    {
        case ElectronScattering::None:
            _com = new NoScatteringHelper(this);
            _numElec = 0;
            break;
        case ElectronScattering::Free:
            _com = new FreeComptonHelper(this);
            _numElec = numAtoms;
            break;
        case ElectronScattering::FreeWithPolarization:
            _com = new FreeComptonWithPolarizationHelper(this);
            _numElec = numAtoms;
            break;
    }

    auto radGrid = find<Configuration>()->radiationFieldWLG();
    if (!radGrid->isContiguous()) throw FATALERROR("Radiation field must consist of consecutive wavelength bins");

    // load the Cloudy wavelength grid
    TextInFile optGrid(this, "XRayCloudyGasMix_wav.dat", "Optical wavelength grid", true);
    optGrid.addColumn("Cloudy wavelength grid", "wavelength", "Ry");
    Array borders = optGrid.readAllColumns()[0];
    vector<double> borderv(std::begin(borders), std::end(borders));
    _opticalWavelengthGrid = new ListBorderWavelengthGrid(this, borderv, true, true);

    // load the emission lines
    TextInFile linesFile(this, "XRayCloudyGasMix_lines.dat", "Cloudy lines", true);
    linesFile.addColumn("mass", "mass", "amu");
    linesFile.addColumn("center", "wavelength", "Angstrom");
    auto lines = linesFile.readAllColumns();
    _lineMassv = lines[0];
    _lineCenterv = lines[1];

    // configure Cloudy
    _cloudyConfig.rad.numBins = radGrid->numBins();
    _cloudyConfig.rad.borderv = Constants::iRy() / radGrid->borderv();  // m -> Ry
    _cloudyConfig.rad.widthv = radGrid->dlambdav();
    _cloudyConfig.rad.minFlux = radMin();
    _cloudyConfig.wav.numBins = _opticalWavelengthGrid->numBins();
    _cloudyConfig.wav.lambdav = _opticalWavelengthGrid->lambdav();
    _cloudyConfig.numIons = numIons;
    _cloudyConfig.numLines = _lineCenterv.size();
    _cloudyConfig.ionIndex = [](const string& name) {
        auto ion = AtomUtils::parseIon(name);
        return AtomUtils::ionIndex(ion.first, ion.second);
    };

    // a plain name is looked up on the PATH by the shell; any other path must point to an existing file and is made
    // absolute, because Cloudy is invoked from within its run directory rather than from the current directory
    string execPath = cloudyExecPath();
    bool isPlainName = !StringUtils::contains(execPath, "/") && !StringUtils::contains(execPath, "\\");
    if (!isPlainName)
    {
        if (!System::isFile(execPath)) throw FATALERROR("The Cloudy executable does not exist: " + execPath);
        if (!StringUtils::isAbsolutePath(execPath)) execPath = System::canonicalPath(execPath);
        if (!StringUtils::isAbsolutePath(execPath))
            throw FATALERROR("Could not determine the absolute path of the Cloudy executable: " + cloudyExecPath());
    }
    _cloudyConfig.execPath = execPath;

    CloudyWrapper::LookupOptions lookup;
    switch (tableMode())
    {
        case TableMode::Direct: lookup.mode = CloudyWrapper::LookupOptions::Mode::Direct; break;
        case TableMode::Record: lookup.mode = CloudyWrapper::LookupOptions::Mode::Record; break;
        case TableMode::Lookup: lookup.mode = CloudyWrapper::LookupOptions::Mode::Lookup; break;
    }
    lookup.numNeighbors = numNeighbors();
    lookup.matchDistance = matchDistance();
    lookup.interpolationDistance = interpolationDistance();
    if (tableMode() == TableMode::Lookup && interpolationDistance() < matchDistance())
        throw FATALERROR("The interpolation distance must not be smaller than the match distance");

    // a relative table directory is relative to the output directory, independent of the current directory
    _tablePath = tableDirectory();
    if (!StringUtils::isAbsolutePath(_tablePath))
        _tablePath = StringUtils::joinPaths(find<FilePaths>()->outputPath(), _tablePath);

    // each process uses its own run directories
    if (!System::makeDir(_tablePath)) throw FATALERROR("Could not create the Cloudy table directory: " + _tablePath);
    string runsName = ProcessManager::isMultiProc() ? "runs_" + StringUtils::toString(ProcessManager::rank()) : "runs";
    string runsPath = StringUtils::joinPaths(_tablePath, runsName);

    _cloudyWrapper = new CloudyWrapper;
    _cloudyWrapper->setup(_cloudyConfig, FilePaths::resource("XRayCloudyGasMix_template.in"), _tablePath, runsPath,
                          lookup);

    if (_cloudyWrapper->usesTable())
        find<Log>()->info("Loaded Cloudy lookup table with " + std::to_string(_cloudyWrapper->tableSize())
                          + " runs from " + _tablePath);
}

////////////////////////////////////////////////////////////////////

XRayCloudyGasMix::~XRayCloudyGasMix()
{
    delete _com;
    delete _cloudyWrapper;
}

////////////////////////////////////////////////////////////////////

int XRayCloudyGasMix::indexForLambda(double lambda) const
{
    return NR::locateFail(_opticalWavelengthGrid->borderv(), lambda);
}

////////////////////////////////////////////////////////////////////

Array XRayCloudyGasMix::getKappaScaCum(const MaterialState* state, int ell) const
{
    Array kappaScaCum(_numElec + 1);
    kappaScaCum[0] = 0;
    for (int e = 0; e < _numElec - 1; e++)
    {
        int index = _indexKappaScaCum + ell * (_numElec - 1) + e;
        kappaScaCum[e + 1] = state->custom(index);
    }
    kappaScaCum[_numElec] = 1;

    return kappaScaCum;
}

////////////////////////////////////////////////////////////////////

void XRayCloudyGasMix::setKappaScaCum(MaterialState* state, int ell, const Array& kappaScaCum) const
{
    // kappaScaCum[0] = 0 by definition
    for (int e = 0; e < _numElec - 1; e++)
    {
        int index = _indexKappaScaCum + ell * (_numElec - 1) + e;
        state->setCustom(index, kappaScaCum[e + 1]);
    }
    // kappaScaCum[_numElec] = 1 by definition
}

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType XRayCloudyGasMix::materialType() const
{
    return MaterialMix::MaterialType::Gas;
}

////////////////////////////////////////////////////////////////////

bool XRayCloudyGasMix::hasPolarizedScattering() const
{
    return electronScattering() == ElectronScattering::FreeWithPolarization;
}

////////////////////////////////////////////////////////////////////

bool XRayCloudyGasMix::hasExtraSpecificState() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

MaterialMix::DynamicStateType XRayCloudyGasMix::hasDynamicMediumState() const
{
    return DynamicStateType::Primary;
}

////////////////////////////////////////////////////////////////////

bool XRayCloudyGasMix::hasScatteringDispersion() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool XRayCloudyGasMix::hasContinuumEmission() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool XRayCloudyGasMix::hasLineEmission() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

#define setAbundance(i, value) setCustom(_indexAbundances + (i), (value))
#define getAbundance(i) custom(_indexAbundances + (i))
#define setVTherm(a, value) setCustom(_indexThermalVelocity + (a), (value))
#define getVTherm(a) custom(_indexThermalVelocity + (a))
#define setKappaAbs(ell, value) setCustom(_indexKappaAbs + (ell), (value))
#define getKappaAbs(ell) custom(_indexKappaAbs + (ell))
#define setKappaSca(ell, value) setCustom(_indexKappaSca + (ell), (value))
#define getKappaSca(ell) custom(_indexKappaSca + (ell))
#define setEmissivity(ell, value) setCustom(_indexEmissivity + (ell), (value))
#define getEmissivity(ell) custom(_indexEmissivity + (ell))
#define setLineEmissivity(ell, value) setCustom(_indexLineEmissivity + (ell), (value))
#define getLineEmissivity(ell) custom(_indexLineEmissivity + (ell))

////////////////////////////////////////////////////////////////////

double XRayCloudyGasMix::kappaAbs(const MaterialState* state, int ell) const
{
    return ell >= 0 ? state->getKappaAbs(ell) : 0.;
}

////////////////////////////////////////////////////////////////////

double XRayCloudyGasMix::kappaSca(const MaterialState* state, int ell) const
{
    return ell >= 0 && _numElec > 0 ? state->getKappaSca(ell) : 0.;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> XRayCloudyGasMix::specificStateVariableInfo() const
{
    vector<StateVariable> result{StateVariable::numberDensity(), StateVariable::temperature(),
                                 StateVariable::metallicity()};

    // To save memory here, we could have some system that allows to only allocate memory for non-empty cells.
    // I.e. have these state variables for only cells with non-zero number density.
    // would have to adjust the MediumState class to support this

    // next available custom variable index
    int index = 0;

    const_cast<XRayCloudyGasMix*>(this)->_indexAbundances = index;
    for (int i = 0; i < numIons; i++)
    {
        const auto& ion = _ionParamv[i];
        string name = AtomUtils::ionName(ion.Z, ion.N);

        result.push_back(StateVariable::custom(index++, name + " abundance", "numbervolumedensity"));
    }

    const_cast<XRayCloudyGasMix*>(this)->_indexThermalVelocity = index;
    for (int a = 0; a < numAtoms; a++)
    {
        string Z = StringUtils::toString(a + 1);
        result.push_back(StateVariable::custom(index++, "thermal velocity for Z=" + Z, "velocity"));
    }

    const_cast<XRayCloudyGasMix*>(this)->_indexKappaAbs = index;
    for (int ell = 0; ell < _cloudyConfig.wav.numBins; ell++)
    {
        string lambda = StringUtils::toString(_cloudyConfig.wav.lambdav[ell]);
        result.push_back(StateVariable::custom(index++, "absorption opacity at lambda=" + lambda + " m", "opacity"));
    }

    if (_numElec > 0)
    {
        const_cast<XRayCloudyGasMix*>(this)->_indexKappaSca = index;
        for (int ell = 0; ell < _cloudyConfig.wav.numBins; ell++)
        {
            string lambda = StringUtils::toString(_cloudyConfig.wav.lambdav[ell]);
            result.push_back(
                StateVariable::custom(index++, "scattering opacity at lambda=" + lambda + " m", "opacity"));
        }

        const_cast<XRayCloudyGasMix*>(this)->_indexKappaScaCum = index;
        for (int ell = 0; ell < _cloudyConfig.wav.numBins; ell++)
            for (int e = 0; e < _numElec - 1; e++)
                result.push_back(StateVariable::custom(index++, "cumulative scattering probability", "1"));
    }

    const_cast<XRayCloudyGasMix*>(this)->_indexEmissivity = index;
    for (int ell = 0; ell < _cloudyConfig.wav.numBins; ell++)
        result.push_back(StateVariable::custom(index++, "volume emissivity", "powervolumedensity"));

    const_cast<XRayCloudyGasMix*>(this)->_indexLineEmissivity = index;
    for (int l = 0; l < _cloudyConfig.numLines; l++)
        result.push_back(StateVariable::custom(index++, "line emissivity", "bolluminosityvolumedensity"));

    return result;
}

////////////////////////////////////////////////////////////////////

void XRayCloudyGasMix::initializeSpecificState(MaterialState* state, double metallicity, double /*temperature*/,
                                               const Array& /*params*/) const
{
    state->setMetallicity(metallicity >= 0. ? metallicity : defaultMetallicity());

    // initialize all cells to empty
    storeOutput(state, _cloudyWrapper->empty());
}

////////////////////////////////////////////////////////////////////

UpdateStatus XRayCloudyGasMix::updateSpecificState(MaterialState* state, const Array& Jv) const
{
    Cloudy::Input input;
    input.hden = state->numberDensity();
    input.metal = state->metallicity();
    input.radv = Jv;

    storeOutput(state, _cloudyWrapper->query(input));

    UpdateStatus status;
    status.updateNotConverged();
    return status;
}

////////////////////////////////////////////////////////////////////

bool XRayCloudyGasMix::isSpecificStateConverged(int numCells, int numUpdated, int numNotConverged,
                                                MaterialState* /*currentAggregate*/,
                                                MaterialState* /*previousAggregate*/) const
{
    // log how the cells were updated, summed over all processes
    auto statistics = _cloudyWrapper->takeStatistics();
    Array counts = {static_cast<double>(statistics.numRuns), static_cast<double>(statistics.numMatches),
                    static_cast<double>(statistics.numInterpolations)};
    ProcessManager::sumToAll(counts);
    auto log = find<Log>();
    string message = "  Cloudy runs: " + std::to_string(static_cast<size_t>(counts[0]));
    if (tableMode() == TableMode::Lookup)
    {
        message += ", reused: " + std::to_string(static_cast<size_t>(counts[1]))
                   + ", interpolated: " + std::to_string(static_cast<size_t>(counts[2]));
    }
    log->info(message);

    if (_cloudyWrapper->usesTable())
    {
        // share the new runs between processes, so that all processes hold all runs
        if (ProcessManager::isMultiProc())
        {
            ProcessManager::broadcastAllToAll([this](vector<double>& data) { _cloudyWrapper->exportNewRuns(data); },
                                              [this](const vector<double>& data) { _cloudyWrapper->importRuns(data); });
        }

        // save the lookup table; all queries have completed at this point
        if (ProcessManager::isRoot() && _cloudyWrapper->save())
        {
            log->info("  Saved Cloudy lookup table with " + std::to_string(_cloudyWrapper->tableSize()) + " runs to "
                      + _tablePath);
        }
    }

    return numCells == numUpdated && numNotConverged == 0;
}

////////////////////////////////////////////////////////////////////

double XRayCloudyGasMix::mass() const
{
    return Constants::Mproton();
}

////////////////////////////////////////////////////////////////////

double XRayCloudyGasMix::sectionAbs(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double XRayCloudyGasMix::sectionSca(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double XRayCloudyGasMix::sectionExt(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double XRayCloudyGasMix::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    return kappaAbs(state, indexForLambda(lambda));
}

////////////////////////////////////////////////////////////////////

double XRayCloudyGasMix::opacitySca(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    return kappaSca(state, indexForLambda(lambda));
}

////////////////////////////////////////////////////////////////////

double XRayCloudyGasMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    int ell = indexForLambda(lambda);
    return kappaAbs(state, ell) + kappaSca(state, ell);
}

////////////////////////////////////////////////////////////////////

void XRayCloudyGasMix::setScatteringInfoIfNeeded(PhotonPacket::ScatteringInfo* scatinfo, double lambda,
                                                 const MaterialState* state) const
{
    if (!scatinfo->valid)
    {
        scatinfo->valid = true;

        // scattering only happens inside the wavelength grid, but clip to be safe
        int ell = NR::locateClip(_opticalWavelengthGrid->borderv(), lambda);
        Array kappaScaCum = getKappaScaCum(state, ell);

        scatinfo->species = NR::locateClip(kappaScaCum, random()->uniform());
        int a = scatinfo->species;
        scatinfo->velocity = state->getVTherm(a) * random()->maxwell();
    }
}

////////////////////////////////////////////////////////////////////

bool XRayCloudyGasMix::peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs,
                                         Direction bfky, const MaterialState* state, const PhotonPacket* pp) const
{
    // draw a random scattering channel and atom velocity, unless a previous peel-off stored this already
    auto scatinfo = const_cast<PhotonPacket*>(pp)->getScatteringInfo();
    setScatteringInfoIfNeeded(scatinfo, lambda, state);

    // only free electron scattering is supported
    int Z = scatinfo->species + 1;

    // scattering in electron rest frame
    lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);
    _com->peeloffScattering(I, Q, U, V, lambda, Z, Z, pp->direction(), bfkobs, bfky, pp);
    lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfkobs, scatinfo->velocity);

    return false;
}

////////////////////////////////////////////////////////////////////

void XRayCloudyGasMix::performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const
{
    // draw a random scattering channel and atom velocity, unless a previous peel-off stored this already
    auto scatinfo = pp->getScatteringInfo();
    setScatteringInfoIfNeeded(scatinfo, lambda, state);

    // only free electron scattering is supported
    int Z = scatinfo->species + 1;

    // scattering in electron rest frame
    lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);
    Direction bfknew = _com->performScattering(lambda, Z, Z, pp->direction(), pp);
    lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfknew, scatinfo->velocity);

    // execute the scattering event in the photon packet
    pp->scatter(bfknew, state->bulkVelocity(), lambda);
}

////////////////////////////////////////////////////////////////////

DisjointWavelengthGrid* XRayCloudyGasMix::emissionWavelengthGrid() const
{
    return _opticalWavelengthGrid;
}

////////////////////////////////////////////////////////////////////

Array XRayCloudyGasMix::emissionSpectrum(const MaterialState* state, const Array& /*Jv*/) const
{
    int numBins = _cloudyConfig.wav.numBins;
    Array emis(numBins + 2);  // zero at both ends
    for (int ell = 0; ell < numBins; ell++) emis[ell + 1] = state->getEmissivity(ell);
    return emis * state->volume();
}

////////////////////////////////////////////////////////////////////

Array XRayCloudyGasMix::lineEmissionCenters() const
{
    return _lineCenterv;
}

////////////////////////////////////////////////////////////////////

Array XRayCloudyGasMix::lineEmissionMasses() const
{
    return _lineMassv;
}

////////////////////////////////////////////////////////////////////

Array XRayCloudyGasMix::lineEmissionSpectrum(const MaterialState* state, const Array& /*Jv*/) const
{
    Array luminosities(_cloudyConfig.numLines);
    for (int l = 0; l != _cloudyConfig.numLines; l++) luminosities[l] = state->getLineEmissivity(l);
    return luminosities * state->volume();
}

////////////////////////////////////////////////////////////////////

double XRayCloudyGasMix::indicativeTemperature(const MaterialState* state, const Array& /*Jv*/) const
{
    return state->temperature();
}

////////////////////////////////////////////////////////////////////

void XRayCloudyGasMix::storeOutput(MaterialState* state, const Cloudy::Output& output) const
{
    const auto& lambdav = _opticalWavelengthGrid->lambdav();

    // temperature
    double temp = output.temp;
    state->setTemperature(temp);

    // thermal velocity
    for (int a = 0; a < numAtoms; a++)
    {
        int Z = a + 1;
        state->setVTherm(a, vtherm(temp, AtomUtils::mass(Z)));
    }

    // abundances
    for (int i = 0; i < numIons; i++) state->setAbundance(i, output.abunv[i]);

    // optical properties
    for (int ell = 0; ell < _cloudyConfig.wav.numBins; ell++)
    {
        // absorption and emission
        state->setKappaAbs(ell, output.opacv[ell]);
        state->setEmissivity(ell, output.emisv[ell]);

        if (_numElec > 0)
        {
            // recalculate scattering with new abundances; only free electron scattering is supported,
            // so the scattering species is the atom and all ions of the same atom are accumulated
            Array kappaScaFractions(0., _numElec);
            for (int i = 0; i < numIons; i++)
            {
                const auto& ion = _ionParamv[i];
                kappaScaFractions[ion.Z - 1] += _com->sectionSca(lambdav[ell], ion.Z, ion.N) * output.abunv[i];
            }

            // determine the normalized cumulative probability distribution and the opacity;
            // without scattering (e.g. an empty cell) the distribution is irrelevant but must remain finite
            double kappaSca = kappaScaFractions.sum();
            Array kappaScaCum;
            if (kappaSca > 0.)
                NR::cdf(kappaScaCum, kappaScaFractions);
            else
                NR::cdf(kappaScaCum, Array(1., _numElec));

            state->setKappaSca(ell, kappaSca);
            setKappaScaCum(state, ell, kappaScaCum);
        }
    }

    // emission lines
    for (int l = 0; l < _cloudyConfig.numLines; l++) state->setLineEmissivity(l, output.linev[l]);
}

////////////////////////////////////////////////////////////////////
