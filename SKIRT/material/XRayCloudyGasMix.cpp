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
#include <unordered_set>

////////////////////////////////////////////////////////////////////

namespace
{
    constexpr int numAtoms = 30;
    constexpr int numIons = 495;

    using IonicCrossSections::nonZeroRange;
    using IonicCrossSections::wavelengthToFromEnergy;

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

    // calculate the cross sections on the Cloudy wavelength grid
    setupCrossSections();

    // configure Cloudy
    _cloudyConfig.rad.numBins = radGrid->numBins();
    _cloudyConfig.rad.borderv = Constants::iRy() / radGrid->borderv();  // m -> Ry
    _cloudyConfig.rad.widthv = radGrid->dlambdav();
    _cloudyConfig.rad.minFlux = radMin();
    _cloudyConfig.wav.numBins = _opticalWavelengthGrid->numBins();
    _cloudyConfig.wav.lambdav = _opticalWavelengthGrid->lambdav();
    _cloudyConfig.numIons = numIons;
    _cloudyConfig.numLines = 0;
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

void XRayCloudyGasMix::setupCrossSections()
{
    int numBins = _opticalWavelengthGrid->numBins();
    const Array& lambdav = _opticalWavelengthGrid->lambdav();

    // the photon energy (eV) at the characteristic wavelength of each bin
    _energyv.resize(numBins);
    for (int ell = 0; ell < numBins; ell++) _energyv[ell] = wavelengthToFromEnergy(lambdav[ell]);

    // load the photo-absorption and fluorescence data for all ions
    std::unordered_set<int> ionSet;
    for (int i = 0; i < numIons; i++) ionSet.insert(i);
    auto paResources = IonicCrossSections::loadPhotoAbsorption(this, ionSet);
    auto fluoResources = IonicCrossSections::loadFluorescence(this, ionSet);

    // the photo-absorption cross section, summed over the subshells of each ion
    _sigmaAbsvv.resize(numIons);
    _sigmaFluvv.resize(numIons);
    _subshellvv.resize(numIons);
    for (const auto& paRes : paResources)
    {
        int i = AtomUtils::ionIndex(paRes.Z, paRes.N);
        auto& sigmaAbsv = _sigmaAbsvv[i];
        if (sigmaAbsv.size() == 0) sigmaAbsv.resize(numBins);
        for (int ell = 0; ell < numBins; ell++) sigmaAbsv[ell] += paRes.photoAbsorbSection(_energyv[ell]);

        // gather the fluorescence lines that may follow a vacancy in this subshell
        vector<double> yields;
        int firstLine = _fluoLinev.size();
        for (const auto& fluoRes : fluoResources)
        {
            if (fluoRes.Z == paRes.Z && fluoRes.N == paRes.N && fluoRes.n == paRes.n && fluoRes.l == paRes.l
                && fluoRes.omega > 0.)
            {
                _fluoLinev.push_back({fluoRes.Z, wavelengthToFromEnergy(fluoRes.E), fluoRes.E, fluoRes.W / 2.});
                yields.push_back(fluoRes.omega);
            }
        }
        if (yields.empty()) continue;

        // store the subshell with its line distribution, which does not depend on wavelength or cell
        Subshell subshell{paRes, 0., firstLine, Array()};
        subshell.yield = NR::cdf(subshell.cumLinev, Array(yields.data(), yields.size()));
        _subshellvv[i].push_back(subshell);

        // the fluorescence part of the cross section, i.e. the cross section times the total yield
        auto& sigmaFluv = _sigmaFluvv[i];
        if (sigmaFluv.size() == 0) sigmaFluv.resize(numBins);
        for (int ell = 0; ell < numBins; ell++)
            sigmaFluv[ell] += paRes.photoAbsorbSection(_energyv[ell]) * subshell.yield;
    }

    // the Compton cross section per electron
    if (_numElec > 0)
    {
        _sigmaElecv.resize(numBins);
        for (int ell = 0; ell < numBins; ell++) _sigmaElecv[ell] = _com->sectionSca(lambdav[ell], 1, 1);
    }

    find<Log>()->info("Loaded photo-absorption data for " + std::to_string(paResources.size())
                      + " subshells and " + std::to_string(_fluoLinev.size()) + " fluorescence lines");
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

bool XRayCloudyGasMix::scatteringEmulatesSecondaryEmission() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool XRayCloudyGasMix::hasContinuumEmission() const
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
#define setKappaFlu(ell, value) setCustom(_indexKappaFlu + (ell), (value))
#define getKappaFlu(ell) custom(_indexKappaFlu + (ell))
#define setElectronDensity(value) setCustom(_indexElectronDensity, (value))
#define getElectronDensity() custom(_indexElectronDensity)
#define setElementCum(a, value) setCustom(_indexElementCum + (a), (value))
#define getElementCum(a) custom(_indexElementCum + (a))
#define setEmissivity(ell, value) setCustom(_indexEmissivity + (ell), (value))
#define getEmissivity(ell) custom(_indexEmissivity + (ell))

////////////////////////////////////////////////////////////////////

double XRayCloudyGasMix::kappaAbs(const MaterialState* state, int ell) const
{
    return ell >= 0 ? state->getKappaAbs(ell) : 0.;
}

////////////////////////////////////////////////////////////////////

double XRayCloudyGasMix::kappaElec(const MaterialState* state, int ell) const
{
    return ell >= 0 && _numElec > 0 ? state->getElectronDensity() * _sigmaElecv[ell] : 0.;
}

////////////////////////////////////////////////////////////////////

double XRayCloudyGasMix::kappaFlu(const MaterialState* state, int ell) const
{
    return ell >= 0 ? state->getKappaFlu(ell) : 0.;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> XRayCloudyGasMix::specificStateVariableInfo() const
{
    vector<StateVariable> result{StateVariable::numberDensity(), StateVariable::temperature(),
                                 StateVariable::metallicity()};

    // next available custom variable index
    int index = 0;
    auto self = const_cast<XRayCloudyGasMix*>(this);
    int numBins = _cloudyConfig.wav.numBins;

    self->_indexAbundances = index;
    for (int i = 0; i < numIons; i++)
    {
        const auto& ion = _ionParamv[i];
        string name = AtomUtils::ionName(ion.Z, ion.N);

        result.push_back(StateVariable::custom(index++, name + " abundance", "numbervolumedensity"));
    }

    self->_indexThermalVelocity = index;
    for (int a = 0; a < numAtoms; a++)
    {
        string Z = StringUtils::toString(a + 1);
        result.push_back(StateVariable::custom(index++, "thermal velocity for Z=" + Z, "velocity"));
    }

    self->_indexKappaAbs = index;
    for (int ell = 0; ell < numBins; ell++)
    {
        string lambda = StringUtils::toString(_cloudyConfig.wav.lambdav[ell]);
        result.push_back(StateVariable::custom(index++, "absorption opacity at lambda=" + lambda + " m", "opacity"));
    }

    self->_indexKappaFlu = index;
    for (int ell = 0; ell < numBins; ell++)
    {
        string lambda = StringUtils::toString(_cloudyConfig.wav.lambdav[ell]);
        result.push_back(
            StateVariable::custom(index++, "fluorescence opacity at lambda=" + lambda + " m", "opacity"));
    }

    if (_numElec > 0)
    {
        self->_indexElectronDensity = index;
        result.push_back(StateVariable::custom(index++, "density of scattering electrons", "numbervolumedensity"));

        self->_indexElementCum = index;
        for (int a = 0; a < _numElec - 1; a++)
            result.push_back(StateVariable::custom(index++, "cumulative electron scattering probability", ""));
    }

    self->_indexEmissivity = index;
    for (int ell = 0; ell < numBins; ell++)
        result.push_back(StateVariable::custom(index++, "volume emissivity", "wavelengthmonluminosityvolumedensity"));

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
    int ell = indexForLambda(lambda);
    return kappaElec(state, ell) + kappaFlu(state, ell);
}

////////////////////////////////////////////////////////////////////

double XRayCloudyGasMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    int ell = indexForLambda(lambda);
    return kappaAbs(state, ell) + kappaElec(state, ell) + kappaFlu(state, ell);
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

        // choose between electron scattering and fluorescence in proportion to their opacities
        double kappaElecValue = kappaElec(state, ell);
        double kappaFluValue = kappaFlu(state, ell);
        if (kappaFluValue <= 0. || random()->uniform() * (kappaElecValue + kappaFluValue) < kappaElecValue)
        {
            // electron scattering: choose the atom whose electrons scatter, which determines the thermal velocity
            Array elementCum(_numElec + 1);
            for (int a = 1; a < _numElec; a++) elementCum[a] = state->getElementCum(a - 1);
            elementCum[_numElec] = 1.;

            int a = NR::locateClip(elementCum, random()->uniform());
            scatinfo->species = a;
            scatinfo->velocity = state->getVTherm(a) * random()->maxwell();
        }
        else
        {
            // fluorescence: choose the ion in proportion to its fluorescence opacity, falling back to the last
            // contributing ion should round-off errors exhaust the total
            double x = random()->uniform() * kappaFluValue;
            int ion = -1;
            for (int i = 0; i < numIons; i++)
            {
                if (_sigmaFluvv[i].size() == 0) continue;
                double kappa = state->getAbundance(i) * _sigmaFluvv[i][ell];
                if (kappa <= 0.) continue;
                ion = i;
                x -= kappa;
                if (x < 0.) break;
            }
            if (ion < 0) throw FATALERROR("No ion contributes to the fluorescence opacity");

            // choose the subshell in proportion to its cross section times its total yield, as in the fluorescence
            // cross section, and then the line from the distribution of the subshell
            const auto& subshells = _subshellvv[ion];
            double E = _energyv[ell];
            Array weights(subshells.size());
            for (size_t s = 0; s != subshells.size(); ++s)
                weights[s] = subshells[s].pa.photoAbsorbSection(E) * subshells[s].yield;
            Array cumWeights;
            NR::cdf(cumWeights, weights);
            const auto& subshell = subshells[NR::locateClip(cumWeights, random()->uniform())];
            int f = subshell.firstLine + NR::locateClip(subshell.cumLinev, random()->uniform());
            const auto& line = _fluoLinev[f];

            scatinfo->species = numAtoms + f;
            scatinfo->velocity = state->getVTherm(line.Z - 1) * random()->maxwell();

            if (line.width == 0.)
            {
                // for a zero-width line, simply copy the central wavelength
                scatinfo->lambda = line.lambda;
            }
            else
            {
                // otherwise sample a wavelength from the Lorentz line shape in energy space; the tails of the Lorentz
                // distribution are very long, occasionally resulting in negative energies; therefore we loop until
                // the sampled wavelength is meaningful
                while (true)
                {
                    scatinfo->lambda = wavelengthToFromEnergy(line.energy + line.width * random()->lorentz());
                    if (nonZeroRange.contains(scatinfo->lambda)) break;
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

bool XRayCloudyGasMix::peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs,
                                         Direction bfky, const MaterialState* state, const PhotonPacket* pp) const
{
    // draw a random scattering channel and atom velocity, unless a previous peel-off stored this already
    auto scatinfo = const_cast<PhotonPacket*>(pp)->getScatteringInfo();
    setScatteringInfoIfNeeded(scatinfo, lambda, state);

    // electron scattering in the electron rest frame
    if (scatinfo->species < numAtoms)
    {
        int Z = scatinfo->species + 1;
        lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);
        _com->peeloffScattering(I, Q, U, V, lambda, Z, Z, pp->direction(), bfkobs, bfky, pp);
        lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfkobs, scatinfo->velocity);
        return false;
    }

    // fluorescence: unpolarized isotropic emission; the bias weight is trivially 1 and there is no contribution to
    // Q, U, V
    I = 1.;
    lambda = PhotonPacket::shiftedEmissionWavelength(scatinfo->lambda, bfkobs, scatinfo->velocity);
    return true;
}

////////////////////////////////////////////////////////////////////

void XRayCloudyGasMix::performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const
{
    // draw a random scattering channel and atom velocity, unless a previous peel-off stored this already
    auto scatinfo = pp->getScatteringInfo();
    setScatteringInfoIfNeeded(scatinfo, lambda, state);

    Direction bfknew;

    // electron scattering in the electron rest frame
    if (scatinfo->species < numAtoms)
    {
        int Z = scatinfo->species + 1;
        lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);
        bfknew = _com->performScattering(lambda, Z, Z, pp->direction(), pp);
        lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfknew, scatinfo->velocity);
    }

    // fluorescence: unpolarized isotropic emission at the wavelength of the line, emulating secondary emission
    else
    {
        pp->setUnpolarized();
        bfknew = random()->direction();
        lambda = PhotonPacket::shiftedEmissionWavelength(scatinfo->lambda, bfknew, scatinfo->velocity);
        pp->setEmulatedSecondaryOrigin(state->mediumIndex());
    }

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

double XRayCloudyGasMix::indicativeTemperature(const MaterialState* state, const Array& /*Jv*/) const
{
    return state->temperature();
}

////////////////////////////////////////////////////////////////////

void XRayCloudyGasMix::storeOutput(MaterialState* state, const Cloudy::Output& output) const
{
    int numBins = _cloudyConfig.wav.numBins;

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

    // photo-absorption and fluorescence opacities from the ion densities
    Array kappaPhotov(numBins);
    Array kappaFluv(numBins);
    for (int i = 0; i < numIons; i++)
    {
        double n = output.abunv[i];
        if (n <= 0.) continue;
        const Array& sigmaAbsv = _sigmaAbsvv[i];
        if (sigmaAbsv.size() != 0)
            for (int ell = 0; ell < numBins; ell++) kappaPhotov[ell] += n * sigmaAbsv[ell];
        const Array& sigmaFluv = _sigmaFluvv[i];
        if (sigmaFluv.size() != 0)
            for (int ell = 0; ell < numBins; ell++) kappaFluv[ell] += n * sigmaFluv[ell];
    }

    // the absorption opacity excludes the fluorescence part, which is treated as scattering; it can be negative where
    // the total fluorescence yield exceeds unity, so that the albedo exceeds unity (as in XRayIonicGasMix)
    for (int ell = 0; ell < numBins; ell++)
    {
        state->setKappaAbs(ell, kappaPhotov[ell] - kappaFluv[ell]);
        state->setKappaFlu(ell, kappaFluv[ell]);
        state->setEmissivity(ell, output.emisv[ell]);
    }

    // electron scattering: with free-electron scattering for all electrons, the share of each atom is independent of
    // wavelength, so the density of the scattering electrons and the cumulative distribution over the atoms suffice
    if (_numElec > 0)
    {
        Array densities(_numElec);
        for (int i = 0; i < numIons; i++)
        {
            int Z = _ionParamv[i].Z;
            densities[Z - 1] += Z * output.abunv[i];
        }

        // without electrons (e.g. an empty cell) the distribution is irrelevant but must remain finite
        Array elementCum;
        double density = densities.sum();
        if (density > 0.)
            NR::cdf(elementCum, densities);
        else
            NR::cdf(elementCum, Array(1., _numElec));

        state->setElectronDensity(density);
        for (int a = 1; a < _numElec; a++) state->setElementCum(a - 1, elementCum[a]);
    }
}

////////////////////////////////////////////////////////////////////
