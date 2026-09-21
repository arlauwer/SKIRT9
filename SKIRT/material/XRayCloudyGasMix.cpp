/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "XRayCloudyGasMix.hpp"
#include "AtomUtils.hpp"
#include "Cloudy.hpp"
#include "CloudyWrapper.hpp"
#include "ComptonPhaseFunction.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "DipolePhaseFunction.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "ListBorderWavelengthGrid.hpp"
#include "MaterialMix.hpp"
#include "MaterialState.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"
#include "UpdateStatus.hpp"
#include <array>

////////////////////////////////////////////////////////////////////

namespace
{
    constexpr int numAtoms = 30;
    constexpr int numIons = 495;

    constexpr double vtherm(double T, double amu)
    {
        return sqrt(Constants::k() / Constants::amu() * T / amu);
    }

    // convert photon energy in Ry to and from wavelength in m (same conversion in both directions)
    template<typename T> constexpr T wavelengthToFromRydberg(T x)
    {
        constexpr double front = Constants::iRy();
        return front / x;
    }
}

////////////////////////////////////////////////////////////////////

void XRayCloudyGasMix::setupSelfBefore()
{
    MaterialMix::setupSelfBefore();

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

    if (!radGrid->isAdjacent()) throw FATALERROR("Radiation field must consist of consecutive wavelength bins");

    // load optical wavelength grid
    TextInFile optGrid(this, "XRayCloudyGasMix_wav.dat", "Optical wavelength grid", true);
    optGrid.addColumn("Cloudy wavelength grid", "wavelength", "Ry");
    Array borders = optGrid.readAllColumns()[0];

    // load emission lines
    TextInFile linesFile(this, "XRayCloudyGasMix_lines.dat", "Cloudy lines", true);
    linesFile.addColumn("mass", "mass", "amu");
    linesFile.addColumn("center", "wavelength", "Angstrom");
    auto lines = linesFile.readAllColumns();

    // --- Wavelength grid ---
    vector<double> borderv(std::begin(borders), std::end(borders));
    _opticalWavelengthGrid = new ListBorderWavelengthGrid(this, borderv, true, true);

    // --- Lines ---
    int numLines = lines[0].size();
    _lineMassv = lines[0];
    _lineCenterv = lines[1];

    // --- Cloudy ---
    _cloudyConfig.setup(*radGrid, *_opticalWavelengthGrid, radMin(), numIons, numLines, cloudyExecPath());
    _cloudyWrapper.setup(&_cloudyConfig, _tableDirectory);
}

////////////////////////////////////////////////////////////////////

XRayCloudyGasMix::~XRayCloudyGasMix()
{
    delete _com;
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
#define setKappaScaCum(ell, index, value) setCustom(_indexKappaScaCum + (ell) * _numElec + (index), (value))
#define getKappaScaCum(ell, index) custom(_indexKappaScaCum + (ell) * _numElec + (index))
#define setEmissivity(ell, value) setCustom(_indexEmissivity + (ell), (value))
#define getEmissivity(ell) custom(_indexEmissivity + (ell))
#define setLineEmissivity(ell, value) setCustom(_indexLineEmissivity + (ell), (value))
#define getLineEmissivity(ell) custom(_indexLineEmissivity + (ell))

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

        result.push_back(StateVariable::custom(index++, name + " abundance", "dimensionless"));
    }

    const_cast<XRayCloudyGasMix*>(this)->_indexThermalVelocity = index;
    for (int a = 0; a < numAtoms; a++)
    {
        string Z = StringUtils::toString(a + 1);
        result.push_back(StateVariable::custom(index++, "thermal velocity for Z=" + Z, "velocity"));
    }

    const_cast<XRayCloudyGasMix*>(this)->_indexKappaAbs = index;
    for (int ell = 0; ell < _cloudyConfig.wav.numBins; ell++)
        result.push_back(StateVariable::custom(index++, "absorption opacity", "opacity"));

    const_cast<XRayCloudyGasMix*>(this)->_indexKappaSca = index;
    for (int ell = 0; ell < _cloudyConfig.wav.numBins; ell++)
        result.push_back(StateVariable::custom(index++, "scattering opacity", "opacity"));

    const_cast<XRayCloudyGasMix*>(this)->_indexKappaScaCum = index;
    for (int ell = 0; ell < _cloudyConfig.wav.numBins; ell++)
        for (int e = 0; e < _numElec + 1; e++)
            result.push_back(StateVariable::custom(index++, "cumulative scattering probability", "1"));

    const_cast<XRayCloudyGasMix*>(this)->_indexEmissivity = index;
    for (int ell = 0; ell < _cloudyConfig.wav.numBins + 2; ell++)
        result.push_back(StateVariable::custom(index++, "volume emissivity", "powervolumedensity"));

    const_cast<XRayCloudyGasMix*>(this)->_indexLineEmissivity = index;
    for (int ell = 0; ell < _cloudyConfig.numLines; ell++)
        result.push_back(StateVariable::custom(index++, "line emissivity", "bolluminosityvolumedensity"));

    return result;
}

////////////////////////////////////////////////////////////////////

void XRayCloudyGasMix::initializeSpecificState(MaterialState* state, double metallicity, double /*temperature*/,
                                               const Array& /*params*/) const
{
    state->setMetallicity(metallicity >= 0. ? metallicity : defaultMetallicity());

    // initialize all cells to empty
    updateSpecificState(state, _cloudyWrapper.empty());
}

////////////////////////////////////////////////////////////////////

UpdateStatus XRayCloudyGasMix::updateSpecificState(MaterialState* state, const Array& J) const
{
    // Array radWidth = config()->radiationFieldWLG()->dlambdav();
    // double ins = (4. * M_PI * J * radWidth).sum();  // W/m2/m/sr -> W/m2 (integrated mean intensity)

    Cloudy::Input input;

    input.hden = state->numberDensity();
    input.metal = state->metallicity();
    input.radv = J;

    // non-const cast
    auto* cloudyWrapper = const_cast<CloudyWrapper*>(&_cloudyWrapper);
    Cloudy::Output output = cloudyWrapper->query(input);

    updateSpecificState(state, output);

    UpdateStatus status;
    status.updateNotConverged();
    return status;
}

////////////////////////////////////////////////////////////////////

bool XRayCloudyGasMix::isSpecificStateConverged(int numCells, int numUpdated, int numNotConverged,
                                                MaterialState* /*currentAggregate*/,
                                                MaterialState* /*previousAggregate*/) const
{
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
    return 0;
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
    int ell = indexForLambda(lambda);
    if (ell < 0 || ell >= _cloudyConfig.wav.numBins)
        return 0.;
    else
        return state->getKappaAbs(ell);
}

////////////////////////////////////////////////////////////////////

double XRayCloudyGasMix::opacitySca(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    int ell = indexForLambda(lambda);
    if (ell < 0 || ell >= _cloudyConfig.wav.numBins)
        return 0.;
    else
        return state->getKappaSca(ell);
}

////////////////////////////////////////////////////////////////////

double XRayCloudyGasMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    int ell = indexForLambda(lambda);
    if (ell < 0 || ell >= _cloudyConfig.wav.numBins)
        return 0.;
    else
        return state->getKappaAbs(ell) + state->getKappaSca(ell);
}

////////////////////////////////////////////////////////////////////

void XRayCloudyGasMix::setScatteringInfoIfNeeded(PhotonPacket::ScatteringInfo* scatinfo, double lambda,
                                                 const MaterialState* state) const
{
    if (!scatinfo->valid)
    {
        scatinfo->valid = true;

        // copy kappaScaCum from the state
        int ell = indexForLambda(lambda);
        Array kappaScaCum(_numElec + 1);
        for (int i = 0; i < _numElec + 1; i++) kappaScaCum[i] = state->getKappaScaCum(ell, i);

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
    // draw a random fluorescence channel and atom velocity, unless a previous peel-off stored this already
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
    Array emis(_cloudyConfig.wav.numBins + 2);  // requires 0 at both ends
    emis[0] = 0.;
    for (int ell = 1; ell < _cloudyConfig.wav.numBins + 1; ell++) emis[ell] = state->getEmissivity(ell);
    emis[_cloudyConfig.wav.numBins + 1] = 0.;
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

void XRayCloudyGasMix::updateSpecificState(MaterialState* state, const Cloudy::Output& output) const
{
    const auto& lamBorderv = _opticalWavelengthGrid->borderv();

    // temperature
    double temp = output.temp;
    state->setTemperature(temp);

    // thermal velocity
    for (int a = 0; a < numAtoms; a++)
    {
        int Z = a + 1;
        double v = vtherm(temp, AtomUtils::mass(Z));
        state->setVTherm(a, v);
    }

    // abundances
    for (int i = 0; i < numIons; i++) state->setAbundance(i, output.abunv[i]);

    // optical properties
    for (int ell = 0; ell < _cloudyConfig.wav.numBins; ell++)
    {
        double lambda = lamBorderv[ell];

        // absorption and emission
        double abs = output.opacv[ell];
        double emi = output.emisv[ell];
        if (std::isnan(abs) || std::isnan(emi)) throw FATALERROR("Cloudy::readOutput found NaN in opac or emis");
        state->setEmissivity(ell, emi);
        state->setKappaAbs(ell, abs);

        // recalculate scattering with new abundances
        Array kappaScaFractions(0., _numElec);
        Array kappaScaCum;

        // for all ions
        for (int i = 0; i < numIons; i++)
        {
            const auto& ion = _ionParamv[i];

            // only free electron scattering is supported
            int species = ion.Z - 1;

            // accumulate abundances of all ions of same Z
            kappaScaFractions[species] += _com->sectionSca(lambda, ion.Z, ion.N) * output.abunv[i];
        }

        // determine the normalized cumulative probability distribution and the cross section
        double kappaSca = NR::cdf(kappaScaCum, kappaScaFractions);

        state->setKappaSca(ell, kappaSca);
        for (int i = 0; i < _numElec + 1; i++)
        {
            state->setKappaScaCum(ell, i, kappaScaCum[i]);
        }
    }

    // emission lines
    for (int l = 0; l < _cloudyConfig.numLines; l++) state->setLineEmissivity(l, output.linev[l]);
}

////////////////////////////////////////////////////////////////////
