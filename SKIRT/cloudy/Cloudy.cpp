#include "Cloudy.hpp"
#include "Atoms.hpp"
#include "Constants.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>

////////////////////////////////////////////////////////////////////

// convert from Cloudy to SKIRT (SI) units
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

    // erg/s -> W
    inline double convertLum(double lum)
    {
        return lum * 1e-7;
    }

    // erg/s/cm2 -> W/m2
    inline void convertInsArray(Array& ins)
    {
        ins *= 1e-3;
    }

    // 1/cm -> 1/m
    inline void convertOpacArray(Array& kappa)
    {
        kappa *= 1e2;
    }

    // Ryd -> m
    inline double convertWavelength(double E)
    {
        return Constants::iRyd() / E;
    }
}

////////////////////////////////////////////////////////////////////

Cloudy::Cloudy(string basePath, const string& inputTemplate, const CloudyConfig& config)
    : _basePath(basePath), _inputTemplate(inputTemplate), _cloudyConfig(config)
{}

////////////////////////////////////////////////////////////////////

void Cloudy::createInput(const Cloudy::Input& input) const
{
    System::makeDir(_basePath);
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
    readLines(output);
}

////////////////////////////////////////////////////////////////////

void Cloudy::createSim(const Input& input) const
{
    double hden = convertDensity(input.hden);
    double metal = convertMetallicity(input.metal);
    double ins = 4. * M_PI * (input.radv * _cloudyConfig.radWidth).sum();
    ins = ins * 1e3;  // W/m2 -> erg/s/cm2

    // create sim.in
    string temp = _inputTemplate;
    temp = StringUtils::replace(temp, "{hden}", StringUtils::toString(hden));
    temp = StringUtils::replace(temp, "{metal}", StringUtils::toString(metal));
    temp = StringUtils::replace(temp, "{ins}", StringUtils::toString(ins));
    std::ofstream sim = System::ofstream(localPath("sim.in"));
    if (!sim.is_open()) throw FATALERROR("Could not open cloudy sim.in file " + _basePath);
    sim.write(temp.data(), temp.size());
    sim.close();
    if (sim.fail()) throw FATALERROR("Could not write to cloudy sim.in file " + _basePath);
}

////////////////////////////////////////////////////////////////////

void Cloudy::createSed(const Input& input) const
{
    // sed.in
    std::ofstream sed = System::ofstream(localPath("sed.in"));
    if (!sed.is_open()) throw FATALERROR("Could not open cloudy sed.in file " + _basePath);
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
    if (sed.fail()) throw FATALERROR("Could not write to cloudy sed.in file " + _basePath);
}

////////////////////////////////////////////////////////////////////

void Cloudy::readTemp(Output& output) const
{
    string line;

    // temperature
    std::ifstream ovr = System::ifstream(localPath("sim.ovr"));
    if (!ovr.is_open()) throw FATALERROR("Could not open cloudy sim.ovr file " + _basePath);
    getline(ovr, line);  // skip header
    getline(ovr, line);
    ovr.close();
    output.temp = StringUtils::toDouble(StringUtils::split(line, "\t")[1]);
}

////////////////////////////////////////////////////////////////////

void Cloudy::readAbun(Output& output) const
{
    output.abunv.resize(_cloudyConfig.numIons);

    string header;
    string line;

    // abundances
    std::ifstream file = System::ifstream(localPath("sim.species"));
    if (!file.is_open()) throw FATALERROR("Could not open cloudy sim.species file " + _basePath);
    getline(file, header);
    getline(file, line);
    file.close();
    if (file.fail()) throw FATALERROR("Could not read from cloudy sim.species file " + _basePath);

    auto speciesHeader = StringUtils::split(header, "\t");
    auto speciesData = StringUtils::split(line, "\t");
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

void Cloudy::binSegments(const string& fileName, Array& binnedSpectrum, size_t col) const
{
    vector<double> segLambdav, segSpectrumv;  // at columns (0, col)

    string line;

    std::ifstream file = System::ifstream(localPath(fileName));
    if (!file.is_open()) throw FATALERROR("Could not open cloudy " + fileName + " file");
    getline(file, line);  // skip header
    while (getline(file, line))
    {
        auto cols = StringUtils::split(line, "\t");
        if (cols.size() <= col) continue;

        segLambdav.push_back(StringUtils::toDouble(cols[0]));
        segSpectrumv.push_back(StringUtils::toDouble(cols[col]));
    }
    file.close();
    int numFileLambda = segLambdav.size();

    // convert x-y segments into binned spectrum

    // reversed order since Cloudy uses ascending energy
    int numLambdaBins = _cloudyConfig.numLambdaBins;
    binnedSpectrum.resize(numLambdaBins, 0.);

    int segIndex = 0;
    int binIndex = 0;

    while (segIndex < numFileLambda - 1 && binIndex < numLambdaBins)
    {
        int segRevIndex = numFileLambda - 1 - segIndex;
        double x1 = segLambdav[segRevIndex], x2 = segLambdav[segRevIndex - 1];
        double y1 = segSpectrumv[segRevIndex], y2 = segSpectrumv[segRevIndex - 1];
        double b1 = _cloudyConfig.lambdaBorderv[binIndex], b2 = _cloudyConfig.lambdaBorderv[binIndex + 1];

        // convert wavelength to m
        x1 = convertWavelength(x1);
        x2 = convertWavelength(x2);

        // skip non-overlapping cases
        if (x2 <= b1)
        {
            segIndex++;
            continue;
        }
        if (b2 <= x1)
        {
            binIndex++;
            continue;
        }

        // compute overlap and accumulate
        double a = max(x1, b1), b = min(x2, b2);
        double area = std::sqrt(y1 * y2) * (b - a);
        binnedSpectrum[binIndex] += area;

        // advance whichever ends first
        if (x2 < b2)
            segIndex++;  // segment ends first
        else if (b2 < x2)
            binIndex++;  // bin ends first
        else
        {
            segIndex++;
            binIndex++;
        }  // both end at same point
    }
}

////////////////////////////////////////////////////////////////////

void Cloudy::readOpac(Output& output) const
{
    binSegments("sim.opac", output.opacv, 2);
    convertOpacArray(output.opacv);
}

////////////////////////////////////////////////////////////////////

void Cloudy::readEmis(Output& output) const
{
    binSegments("sim.emis", output.emisv, 3);
    convertInsArray(output.emisv);
}

////////////////////////////////////////////////////////////////////

void Cloudy::readLines(Output& output) const
{
    int numLines = _cloudyConfig.numLines;
    output.linev.resize(numLines);

    string line;

    std::ifstream lines = System::ifstream(localPath("sim.lines"));
    getline(lines, line);  // skip header
    getline(lines, line);  // single row

    auto cols = StringUtils::split(line, "\t");

    if ((int)cols.size() != numLines + 1) throw FATALERROR("Cloudy::readLines read wrong number of values");

    // skirt first col (depth)
    for (int i = 1; i < _cloudyConfig.numLines; i++)
    {
        output.linev[i] = convertLum(StringUtils::toDouble(cols[i]));
    }
}

////////////////////////////////////////////////////////////////////

string Cloudy::localPath(const string& filename) const
{
    return StringUtils::joinPaths(_basePath, filename);
}

////////////////////////////////////////////////////////////////////
