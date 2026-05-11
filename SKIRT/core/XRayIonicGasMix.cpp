/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "XRayIonicGasMix.hpp"
#include "Atoms.hpp"
#include "Cloudy.hpp"
#include "CloudyWrapper.hpp"
#include "ComptonPhaseFunction.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "DipolePhaseFunction.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "FilePaths.hpp"
#include "ListBorderWavelengthGrid.hpp"
#include "MaterialMix.hpp"
#include "MaterialState.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "TextInFile.hpp"
#include "UpdateStatus.hpp"
#include <array>

////////////////////////////////////////////////////////////////////

namespace
{
    constexpr int numAtoms = CloudyConfig::numAtoms;
    constexpr int numIons = CloudyConfig::numIons;
    constexpr int numInter = CloudyConfig::numIons;  // Compton scattering
    const auto ions = Atoms::initIons<numIons, numAtoms>();

    constexpr double vtherm(double T, double amu)
    {
        return sqrt(Constants::k() / Constants::amu() * T / amu);
    }

    // convert photon energy in eV to and from wavelength in m (same conversion in both directions)
    constexpr double wavelengthToFromEnergy(double x)
    {
        constexpr double front = Constants::h() * Constants::c() / Constants::Qelectron();
        return front / x;
    }

    // convert photon energy in Ryd to and from wavelength in m (same conversion in both directions)
    template<typename T> constexpr T wavelengthToFromRydberg(T x)
    {
        constexpr double front = Constants::iRyd();
        return front / x;
    }
}

////////////////////////////////////////////////////////////////////

// ---- base class for scattering helpers ---- //

class XRayIonicGasMix::ScatteringHelper
{
public:
    virtual ~ScatteringHelper() {}

    // return scattering cross section for atom in m2
    virtual double sectionSca(double lambda, int Z) const = 0;

    // peel-off unpolarized scattering event: override this in helpers that don't support polarization
    virtual void peeloffScattering(double& /*I*/, double& /*lambda*/, int /*Z*/, Direction /*bfk*/,
                                   Direction /*bfkobs*/) const
    {
        // default implementation does nothing
    }

    // perform unpolarized scattering event: override this in helpers that don't support polarization
    virtual Direction performScattering(double& /*lambda*/, int /*Z*/, Direction /*bfk*/) const
    {
        // default implementation returns null vector
        return Direction();
    }

    // peel-off polarized scattering event: override this in helpers that do support polarization
    virtual void peeloffScattering(double& I, double& /*Q*/, double& /*U*/, double& /*V*/, double& lambda, int Z,
                                   Direction bfk, Direction bfkobs, Direction /*bfky*/,
                                   const StokesVector* /*sv*/) const
    {
        // default implementation calls unpolarized version
        peeloffScattering(I, lambda, Z, bfk, bfkobs);
    }

    // perform polarized scattering event: override this in helpers that do support polarization
    virtual Direction performScattering(double& lambda, int Z, Direction bfk, StokesVector* /*sv*/) const
    {
        // default implementation calls unpolarized version
        return performScattering(lambda, Z, bfk);
    }
};

////////////////////////////////////////////////////////////////////

// ---- no scattering helper ---- //

namespace
{
    // this helper does nothing; it is used as a stub in case there is no scattering of a given type
    class NoScatteringHelper : public XRayIonicGasMix::ScatteringHelper
    {
    public:
        NoScatteringHelper(SimulationItem* /*item*/) {}

        double sectionSca(double /*lambda*/, int /*Z*/) const override { return 0.; }
    };
}

////////////////////////////////////////////////////////////////////

// ---- free-electron Compton scattering helper ---- //

namespace
{
    // transition wavelength from Compton to Thomson scattering
    constexpr double comptonWL = wavelengthToFromEnergy(100.);  // 0.1 keV or 12.4 nm

    // this helper forwards all calls to an external helper class for regular Compton scattering
    // (or Thomson scattering for lower energies, because Compton becomes numerically unstable)
    class FreeComptonHelper : public XRayIonicGasMix::ScatteringHelper
    {
    private:
        ComptonPhaseFunction _cpf;
        DipolePhaseFunction _dpf;

    public:
        FreeComptonHelper(SimulationItem* item)
        {
            auto random = item->find<Random>();
            _cpf.initialize(random);
            _dpf.initialize(random);
        }

        double sectionSca(double lambda, int Z) const override
        {
            double sigma = Z * Constants::sigmaThomson();
            if (lambda < comptonWL) sigma *= _cpf.sectionSca(lambda);
            return sigma;
        }

        void peeloffScattering(double& I, double& lambda, int /*Z*/, Direction bfk, Direction bfkobs) const override
        {
            if (lambda < comptonWL)
            {
                double Q, U, V;
                _cpf.peeloffScattering(I, Q, U, V, lambda, bfk, bfkobs, Direction(), nullptr);
            }
            else
            {
                double Q, U, V;
                _dpf.peeloffScattering(I, Q, U, V, bfk, bfkobs, Direction(), nullptr);
            }
        }

        Direction performScattering(double& lambda, int /*Z*/, Direction bfk) const override
        {
            return lambda < comptonWL ? _cpf.performScattering(lambda, bfk, nullptr)
                                      : _dpf.performScattering(bfk, nullptr);
        }
    };
}

////////////////////////////////////////////////////////////////////

// ---- free-electron Compton with polarization scattering helper ---- //

namespace
{
    // this helper forwards all calls to an external helper class for Compton scattering
    // (or Thomson scattering for lower energies) with support for polarization
    class FreeComptonWithPolarizationHelper : public XRayIonicGasMix::ScatteringHelper
    {
    private:
        ComptonPhaseFunction _cpf;
        DipolePhaseFunction _dpf;

    public:
        FreeComptonWithPolarizationHelper(SimulationItem* item)
        {
            auto random = item->find<Random>();
            _cpf.initialize(random, true);
            _dpf.initialize(random, true);
        }

        double sectionSca(double lambda, int Z) const override
        {
            double sigma = Z * Constants::sigmaThomson();
            if (lambda < comptonWL) sigma *= _cpf.sectionSca(lambda);
            return sigma;
        }

        void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, int /*Z*/, Direction bfk,
                               Direction bfkobs, Direction bfky, const StokesVector* sv) const override
        {
            lambda < comptonWL ? _cpf.peeloffScattering(I, Q, U, V, lambda, bfk, bfkobs, bfky, sv)
                               : _dpf.peeloffScattering(I, Q, U, V, bfk, bfkobs, bfky, sv);
        }

        Direction performScattering(double& lambda, int /*Z*/, Direction bfk, StokesVector* sv) const override
        {
            return lambda < comptonWL ? _cpf.performScattering(lambda, bfk, sv) : _dpf.performScattering(bfk, sv);
        }
    };
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMix::setupSelfBefore()
{
    MaterialMix::setupSelfBefore();

    _log = find<Log>();

    // create scattering helpers depending on the user-configured implementation type;
    // the respective helper constructors load the required bound-electron scattering resources
    switch (electronScattering())
    {
        case ElectronScattering::None: _com = new NoScatteringHelper(this); break;
        case ElectronScattering::Free: _com = new FreeComptonHelper(this); break;
        case ElectronScattering::FreeWithPolarization: _com = new FreeComptonWithPolarizationHelper(this); break;
    }

    // have to set this up before we can call setupCloudyConfig()
    // but setupCloudyConfig needs to be in the setupSelfBefore() function
    opticalWavelengthGrid()->setup();

    setupCloudyConfig();

    // continuum emission wavelength grid
    vector<double> borders(std::begin(_cloudyConfig.lambdaBorderv), std::end(_cloudyConfig.lambdaBorderv));
    _emissionWavelengthGrid = new ListBorderWavelengthGrid(this, borders, true, true);

    // cloudy wrapper
    string basePath = StringUtils::dirPath(FilePaths::resource("XRayIonicGasMix_template.in"));
    _cloudyWrapper.setup(_cloudyConfig, basePath);

    // write backup template
    auto filePaths = find<FilePaths>();
    string templateBackupPath = filePaths->output("template.in");
    std::ofstream templateBackup(templateBackupPath);
    templateBackup << _cloudyWrapper.templateContent();
    templateBackup.close();
}

////////////////////////////////////////////////////////////////////

XRayIonicGasMix::~XRayIonicGasMix()
{
    delete _com;
}

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType XRayIonicGasMix::materialType() const
{
    return MaterialMix::MaterialType::Gas;
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::hasPolarizedScattering() const
{
    return electronScattering() == ElectronScattering::FreeWithPolarization;
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::hasExtraSpecificState() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

MaterialMix::DynamicStateType XRayIonicGasMix::hasDynamicMediumState() const
{
    return DynamicStateType::Primary;
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::hasScatteringDispersion() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::hasContinuumEmission() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::hasLineEmission() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

#define setAbundance(ion, value) setCustom(_indexAbundances + (ion), (value))
#define getAbundance(ion) custom(_indexAbundances + (ion))
#define setVTherm(Z, value) setCustom(_indexThermalVelocity + (Z), (value))
#define getVTherm(Z) custom(_indexThermalVelocity + (Z))
#define setKappaAbs(l, value) setCustom(_indexKappaAbs + (l), (value))
#define getKappaAbs(l) custom(_indexKappaAbs + (l))
#define setKappaSca(l, value) setCustom(_indexKappaSca + (l), (value))
#define getKappaSca(l) custom(_indexKappaSca + (l))
#define setKappaScaCum(l, index, value) setCustom(_indexKappaScaCum + (l) * numInter + (index), (value))
#define getKappaScaCum(l, index) custom(_indexKappaScaCum + (l) * numInter + (index))
#define setEmissivity(l, value) setCustom(_indexEmissivity + (l), (value))
#define getEmissivity(l) custom(_indexEmissivity + (l))
#define setLineEmissivity(l, value) setCustom(_indexLineEmissivity + (l), (value))
#define getLineEmissivity(l) custom(_indexLineEmissivity + (l))

////////////////////////////////////////////////////////////////////

vector<StateVariable> XRayIonicGasMix::specificStateVariableInfo() const
{
    vector<StateVariable> result{StateVariable::numberDensity(), StateVariable::temperature(),
                                 StateVariable::metallicity()};

    // To save memory here, we could have some system that allows to only allocate memory for non-empty cells.
    // I.e. have these state variables for only cells with non-zero number density.

    // next available custom variable index
    int index = 0;

    const_cast<XRayIonicGasMix*>(this)->_indexAbundances = index;
    for (int i = 0; i < numIons; i++)
        result.push_back(StateVariable::custom(index++, "abundance", "numbervolumedensity"));

    const_cast<XRayIonicGasMix*>(this)->_indexThermalVelocity = index;
    for (int a = 0; a < numAtoms; a++) result.push_back(StateVariable::custom(index++, "thermal velocity", "velocity"));

    const_cast<XRayIonicGasMix*>(this)->_indexKappaAbs = index;
    for (int l = 0; l < _cloudyConfig.numLambdaBins; l++)
        result.push_back(StateVariable::custom(index++, "absorption opacity", "opacity"));

    const_cast<XRayIonicGasMix*>(this)->_indexKappaSca = index;
    for (int l = 0; l < _cloudyConfig.numLambdaBins; l++)
        result.push_back(StateVariable::custom(index++, "scattering opacity", "opacity"));

    const_cast<XRayIonicGasMix*>(this)->_indexKappaScaCum = index;
    for (int l = 0; l < _cloudyConfig.numLambdaBins; l++)
        for (int i = 0; i < numInter + 1; ++i)
            result.push_back(StateVariable::custom(index++, "cumulative scattering probability", "1"));

    const_cast<XRayIonicGasMix*>(this)->_indexEmissivity = index;
    for (int l = 0; l < _cloudyConfig.numLambdaBins + 2; l++)
        result.push_back(StateVariable::custom(index++, "volume emissivity", "powervolumedensity"));

    const_cast<XRayIonicGasMix*>(this)->_indexLineEmissivity = index;
    for (int l = 0; l < _cloudyConfig.numLines; l++)
        result.push_back(StateVariable::custom(index++, "line emissivity", "bolluminosityvolumedensity"));

    return result;
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMix::initializeSpecificState(MaterialState* state, double metallicity, double /*temperature*/,
                                              const Array& /*params*/) const
{
    state->setMetallicity(metallicity >= 0. ? metallicity : defaultMetallicity());

    // initialize all cells to empty
    updateSpecificState(state, _cloudyWrapper.empty());
}

////////////////////////////////////////////////////////////////////

UpdateStatus XRayIonicGasMix::updateSpecificState(MaterialState* state, const Array& J) const
{
    // Array radWidth = config()->radiationFieldWLG()->dlambdav();
    // double ins = (4. * M_PI * J * radWidth).sum();  // W/m2/m/sr -> W/m2 (integrated mean intensity)

    Cloudy::Input input;

    input.hden = state->numberDensity();
    input.metal = state->metallicity();
    input.radv = J;

    auto* cloudyWrapper = const_cast<CloudyWrapper*>(&_cloudyWrapper);

    Cloudy::Output output = cloudyWrapper->query(input);

    updateSpecificState(state, output);

    UpdateStatus status;
    status.updateNotConverged();
    return status;
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::isSpecificStateConverged(int numCells, int numUpdated, int numNotConverged,
                                               MaterialState* /*currentAggregate*/,
                                               MaterialState* /*previousAggregate*/) const
{
    return numCells == numUpdated && numNotConverged == 0;
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::mass() const
{
    return Constants::Mproton();
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::sectionAbs(double /*lambda*/) const
{
    return 0;
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::sectionSca(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::sectionExt(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    int ell = NR::locateFail(_cloudyConfig.lambdaBorderv, lambda);
    if (ell < 0 || ell >= _cloudyConfig.numLambdaBins)
        return 0.;
    else
        return state->getKappaAbs(ell);
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::opacitySca(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    int ell = NR::locateFail(_cloudyConfig.lambdaBorderv, lambda);
    if (ell < 0 || ell >= _cloudyConfig.numLambdaBins)
        return 0.;
    else
        return state->getKappaSca(ell);
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    int ell = NR::locateFail(_cloudyConfig.lambdaBorderv, lambda);
    if (ell < 0 || ell >= _cloudyConfig.numLambdaBins)
        return 0.;
    else
        return state->getKappaAbs(ell) + state->getKappaSca(ell);
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMix::setScatteringInfoIfNeeded(PhotonPacket::ScatteringInfo* scatinfo, double lambda,
                                                const MaterialState* state) const
{
    if (!scatinfo->valid)
    {
        scatinfo->valid = true;
        int lam = NR::locateClip(_cloudyConfig.lambdaBorderv, lambda);  // this should (almost) never have to clip
        // scattering can only happen if opacity is non-zero, so lambda should be in range of _lambdaC
        // maybe some Doppler shift but a simple clip should be sufficient
        Array kappaScaCum(numInter + 1);
        for (int i = 0; i < numInter + 1; i++) kappaScaCum[i] = state->getKappaScaCum(lam, i);

        scatinfo->species = NR::locateClip(kappaScaCum, random()->uniform());

        int i = scatinfo->species % numIons;
        int Z = ions[i].Z;

        scatinfo->velocity = state->getVTherm(Z) * random()->maxwell();
    }
}

////////////////////////////////////////////////////////////////////

bool XRayIonicGasMix::peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs,
                                        Direction bfky, const MaterialState* state, const PhotonPacket* pp) const
{
    // draw a random scattering channel and atom velocity, unless a previous peel-off stored this already
    auto scatinfo = const_cast<PhotonPacket*>(pp)->getScatteringInfo();
    setScatteringInfoIfNeeded(scatinfo, lambda, state);

    // if we have dispersion, for electron scattering, adjust the incoming wavelength to the electron rest frame
    if (state->temperature() > 0.)
        lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);

    // Compton scattering in electron rest frame; with support for polarization if enabled
    const auto& ion = ions[scatinfo->species];
    _com->peeloffScattering(I, Q, U, V, lambda, ion.Z, pp->direction(), bfkobs, bfky, pp);

    // if we have dispersion, Doppler-shift the outgoing wavelength from the electron rest frame
    if (state->temperature() > 0.) lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfkobs, scatinfo->velocity);

    return false;
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMix::performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const
{
    // draw a random fluorescence channel and atom velocity, unless a previous peel-off stored this already
    auto scatinfo = pp->getScatteringInfo();
    setScatteringInfoIfNeeded(scatinfo, lambda, state);

    // if we have dispersion, for electron scattering, adjust the incoming wavelength to the electron rest frame
    if (state->temperature() > 0.)
        lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);

    // room for the outgoing direction
    Direction bfknew;

    // Compton scattering, with support for polarization if enabled:
    // determine the new propagation direction and wavelength, and if polarized, update the stokes vector
    const auto& ion = ions[scatinfo->species];
    bfknew = _com->performScattering(lambda, ion.Z, pp->direction(), pp);

    // if we have dispersion, Doppler-shift the outgoing wavelength from the electron rest frame
    if (state->temperature() > 0.) lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfknew, scatinfo->velocity);

    // execute the scattering event in the photon packet
    pp->scatter(bfknew, state->bulkVelocity(), lambda);
}

////////////////////////////////////////////////////////////////////

DisjointWavelengthGrid* XRayIonicGasMix::emissionWavelengthGrid() const
{
    return _emissionWavelengthGrid;
}

////////////////////////////////////////////////////////////////////

Array XRayIonicGasMix::emissionSpectrum(const MaterialState* state, const Array& /*Jv*/) const
{
    Array emis(_cloudyConfig.numLambdaBins + 2);  // requires 0 at both ends
    emis[0] = 0.;
    for (int l = 1; l < _cloudyConfig.numLambdaBins + 1; l++) emis[l] = state->getEmissivity(l);
    emis[_cloudyConfig.numLambdaBins + 1] = 0.;
    return emis * state->volume();
}

////////////////////////////////////////////////////////////////////

Array XRayIonicGasMix::lineEmissionCenters() const
{
    return _cloudyConfig.lineEmisCenterv;
}

////////////////////////////////////////////////////////////////////

Array XRayIonicGasMix::lineEmissionMasses() const
{
    return _cloudyConfig.lineMassv;
}

////////////////////////////////////////////////////////////////////

Array XRayIonicGasMix::lineEmissionSpectrum(const MaterialState* state, const Array& /*Jv*/) const
{
    Array luminosities(_cloudyConfig.numLines);
    for (int l = 0; l != _cloudyConfig.numLines; l++) luminosities[l] = state->getLineEmissivity(l);
    return luminosities * state->volume();
}

////////////////////////////////////////////////////////////////////

double XRayIonicGasMix::indicativeTemperature(const MaterialState* state, const Array& /*Jv*/) const
{
    return state->temperature();
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMix::setupCloudyConfig()
{
    auto config = find<Configuration>();
    auto radGrid = config->radiationFieldWLG();

    if (!radGrid->isAdjacent()) throw FATALERROR("Radiation field must consist of consecutive wavelength bins");
    if (!opticalWavelengthGrid()->isAdjacent())
        throw FATALERROR("Optical wavelength grid must consist of consecutive wavelength bins");

    // load optical wavelength grid
    TextInFile linesFile(this, "XRayIonicGasMix_lines.dat", "Cloudy lines", true);
    linesFile.addColumn("mass", "mass", "amu");
    linesFile.addColumn("center", "wavelength", "Angstrom");
    auto lines = linesFile.readAllColumns();

    // --- Radiation field ---
    _cloudyConfig.numBins = radGrid->numBins();
    _cloudyConfig.radEdges = wavelengthToFromRydberg(radGrid->borderv());
    _cloudyConfig.radWidth = radGrid->dlambdav();
    _cloudyConfig.radMin = radMin();

    // --- Optical properties ---
    _cloudyConfig.numLambdaBins = opticalWavelengthGrid()->numBins();
    _cloudyConfig.lambdaBorderv = opticalWavelengthGrid()->borderv();

    // --- Lines ---
    _cloudyConfig.numLines = lines[0].size();
    _cloudyConfig.lineEmisCenterv = lines[1];
    _cloudyConfig.lineMassv = lines[0];

    // --- Cloudy ---
    _cloudyConfig.cloudyExecPath = cloudyExecPath();
    _cloudyConfig.numDims = _cloudyConfig.numBins + 2;  // + 2 for hden and metallicity
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMix::updateSpecificState(MaterialState* state, const Cloudy::Output& output) const
{
    // temperature
    double temp = output.temp;
    state->setTemperature(temp);

    // thermal velocity
    for (int i = 0; i < numAtoms; i++)
    {
        short Z = i + 1;
        double v = vtherm(temp, Atoms::mass(Z));
        state->setVTherm(i, v);
    }

    // abundances
    for (int i = 0; i < numIons; i++) state->setAbundance(i, output.abunv[i]);

    // optical properties
    for (int ell = 0; ell < _cloudyConfig.numLambdaBins; ell++)
    {
        double lambda = _cloudyConfig.lambdaBorderv[ell];
        double abs = output.opacv[ell];
        double emi = output.emisv[ell];

        if (std::isnan(abs) || std::isnan(emi)) throw FATALERROR("Cloudy::readOutput found NaN in opac or emis");

        state->setEmissivity(ell, emi);
        state->setKappaAbs(ell, abs);

        // --- Recalculate scattering with new abundances ---
        // provide temporary array for the non-normalized fluorescence/scattering contributions (at the current wavelength)
        Array kappaScaFractions(numInter);
        Array kappaScaCum;

        // electron scattering
        for (int i = 0; i < numIons; i++)
        {
            double abun = output.abunv[i];
            const auto& ion = ions[i];

            kappaScaFractions[i] = _com->sectionSca(lambda, ion.Z) * abun;
        }

        // determine the normalized cumulative probability distribution and the cross section
        double kappaSca = NR::cdf(kappaScaCum, kappaScaFractions);

        state->setKappaSca(ell, kappaSca);
        for (int i = 0; i < numInter + 1; i++)
        {
            state->setKappaScaCum(ell, i, kappaScaCum[i]);
        }
    }

    // emission lines
    for (int l = 0; l < _cloudyConfig.numLines; l++) state->setLineEmissivity(l, output.linev[l]);
}

////////////////////////////////////////////////////////////////////
