#include "Cloudy.hpp"
#include "Atoms.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <cstddef>

////////////////////////////////////////////////////////////////////

// convert from SKIRT (SI) to Cloudy units
namespace
{

    // 1/cm3 -> 1/m3
    inline double convertDensity(double hden)
    {
        return hden * 1e6;
    }

    // 1 -> 1
    inline double convertMetallicity(double metallicity)
    {
        return metallicity;
    }

    // W/m2 -> erg/s/cm2
    inline double convertIns(double ins)
    {
        return ins * 1e3;
    }

    // 1/cm -> 1/m
    inline double convertOpac(double kappa)
    {
        return kappa * 1e2;
    }
}

////////////////////////////////////////////////////////////////////

Cloudy::Cloudy(string basePath, const string& inputTemplate, const CloudyConfig& config)
    : _basePath(basePath), _inputTemplate(inputTemplate), _cloudyConfig(config)
{}

////////////////////////////////////////////////////////////////////

void Cloudy::createInput(const Cloudy::Input& input) const
{
    createSim(input);
    createSed(input);
}

////////////////////////////////////////////////////////////////////

bool Cloudy::execute()
{
    string cmd = "cd " + _basePath + " && " + _cloudyConfig.cloudyExecPath + " < sim.in > sim.out";
    int status = System::execute(cmd);
    return status == 0;
}

////////////////////////////////////////////////////////////////////

void Cloudy::readOutput(Output& output) const
{
    readTemp(output);
    readAbun(output);
    readOpac(output);
    readEmis(output);
}

////////////////////////////////////////////////////////////////////

void Cloudy::createSim(const Input& input) const
{
    double hden = convertDensity(input.hden);
    double metal = convertMetallicity(input.metal);
    double ins = 4. * M_PI * (input.radv * _cloudyConfig.radWidth).sum();
    ins = convertIns(ins);

    // create sim.in
    string temp = _inputTemplate;
    temp = StringUtils::replace(temp, "{hden}", StringUtils::toString(hden));
    temp = StringUtils::replace(temp, "{metal}", StringUtils::toString(metal));
    temp = StringUtils::replace(temp, "{ins}", StringUtils::toString(ins));
    std::ofstream in = System::ofstream(localPath("sim.in"));
    in.write(temp.data(), temp.size());
    in.close();
}

////////////////////////////////////////////////////////////////////

void Cloudy::createSed(const Input& input) const
{
    // sed.in
    std::ofstream sed = System::ofstream(localPath("sed.in"));
    // SKIRT  has      radv    (4pi W/m2/m)
    // Cloudy requires Jlambda (~   W/m2/m)
    // so no unit conversion is required
    for (size_t i = 0; i < _cloudyConfig.numBins; i++)
    {
        double left = _cloudyConfig.radEdges[i] * 0.9999;
        double right = _cloudyConfig.radEdges[i + 1] * 1.00001;
        double rad = max(input.radv[i], _cloudyConfig.radMin);  // to avoid 0

        sed << left << "\t" << rad << (i == 0 ? " Flambda" : "") << "\n";
        sed << right << "\t" << rad << "\n";
    }
    sed.close();
}

////////////////////////////////////////////////////////////////////

void Cloudy::readTemp(Output& output) const
{
    string header;
    string line;

    // temperature
    std::ifstream ovr = System::ifstream(localPath("sim.ovr"));
    getline(ovr, header);  // skip header
    getline(ovr, line);
    ovr.close();
    output.temp = StringUtils::toDouble(StringUtils::split(line, "\t")[1]);
}

////////////////////////////////////////////////////////////////////

void Cloudy::readAbun(Output& output) const
{
    output.abunv.resize(_cloudyConfig.numIons);

    string header;
    string species;

    // abundances
    std::ifstream file = System::ifstream(localPath("sim.species"));
    getline(file, header);
    getline(file, species);
    file.close();

    auto speciesHeader = StringUtils::split(header, "\t");
    auto speciesData = StringUtils::split(species, "\t");
    for (size_t i = 1; i < speciesData.size() - 5; i++)  // skip columns: 1, -1, -2, -3, -4, -5
    {
        auto speciesName = speciesHeader[i];
        Atoms::Ion ion = Atoms::parseIon(speciesName);
        if (ion.N > 0)  // not fully ionized
        {
            int ionIndex = Atoms::ionIndex(ion);  // linear index
            double abun = convertDensity(StringUtils::toDouble(speciesData[i]));
            output.abunv[ionIndex] = abun;
        }
    }
}

////////////////////////////////////////////////////////////////////

void Cloudy::readOpac(Output& output) const
{
    int numLambda = _cloudyConfig.numLambda;
    output.opacv.resize(numLambda);

    string line;

    // opacities
    std::ifstream opac = System::ifstream(localPath("sim.opac"));
    getline(opac, line);  // skip header
    int i = 0;
    while (getline(opac, line))
    {
        auto cols = StringUtils::split(line, "\t");
        if (cols.size() >= 3)
        {
            if (i > numLambda) throw FATALERROR("Cloudy::readOpac read too many opacities");

            output.opacv[i] = convertOpac(StringUtils::toDouble(cols[2]));
            i++;
        }
    }
    opac.close();

    if (i != numLambda) throw FATALERROR("Cloudy::readOpac read too few opacities");
}

////////////////////////////////////////////////////////////////////

void Cloudy::readEmis(Output& output) const
{
    int numLambda = _cloudyConfig.numLambda;
    output.emisv.resize(numLambda);

    string line;

    // emissivities
    std::ifstream emis = System::ifstream(localPath("sim.con"));
    getline(emis, line);  // skip header
    int i = 0;
    while (getline(emis, line))
    {
        auto cols = StringUtils::split(line, "\t");
        if (cols.size() >= 4)
        {
            if (i > numLambda) throw FATALERROR("Cloudy::readEmis read too many emissivities");

            output.emisv[i] = convertIns(StringUtils::toDouble(cols[3]));
            i++;
        }
    }
    emis.close();

    if (i != numLambda) throw FATALERROR("Cloudy::readEmis read too few emissivities");
}

////////////////////////////////////////////////////////////////////

string Cloudy::localPath(const string& filename) const
{
    return StringUtils::joinPaths(_basePath, filename);
}

////////////////////////////////////////////////////////////////////
