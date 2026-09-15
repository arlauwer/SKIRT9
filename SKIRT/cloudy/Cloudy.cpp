#include "Cloudy.hpp"
#include "AtomUtils.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <fstream>

////////////////////////////////////////////////////////////////////

namespace
{
    // returns the tab separated fields of the given line, or throws if there are too few.
    vector<string> splitRow(const string& line, size_t minCols, const string& what)
    {
        auto cols = StringUtils::split(line, "\t");
        if (cols.size() < minCols)
            throw FATALERROR("Cloudy: expected at least " + StringUtils::toString(static_cast<int>(minCols))
                             + " columns in " + what);
        return cols;
    }

    // opens a Cloudy output file for reading, or throws.
    std::ifstream openInput(const string& path)
    {
        std::ifstream in = System::ifstream(path);
        if (!in.is_open()) throw FATALERROR("Cloudy: could not open " + path);
        return in;
    }

    // opens a Cloudy input file for writing, or throws.
    std::ofstream openOutput(const string& path)
    {
        std::ofstream out = System::ofstream(path);
        if (!out.is_open()) throw FATALERROR("Cloudy: could not open " + path);
        return out;
    }

    // closes a written file and throws if anything went wrong along the way.
    void closeChecked(std::ofstream& out, const string& path)
    {
        out.close();
        if (out.fail()) throw FATALERROR("Cloudy: could not write " + path);
    }
}

////////////////////////////////////////////////////////////////////

void Cloudy::Output::resize(const CloudyConfig& config)
{
    temp = 0.;
    abunv.resize(config.numIons, 0.);
    opacv.resize(config.wav.numBins, 0.);
    emisv.resize(config.wav.numBins, 0.);
    linev.resize(config.numLines, 0.);
}

////////////////////////////////////////////////////////////////////

const vector<int>& CloudySpeciesHeader::columnsFor(const vector<string>& headerCols)
{
    std::unique_lock<std::mutex> lock(_mutex);

    if (_header.empty())
    {
        // first run: resolve every column name once
        // column 0 is the depth, the last five columns are Cloudy bookkeeping
        const size_t numTrailing = 5;
        if (headerCols.size() <= 1 + numTrailing) throw FATALERROR("Cloudy: the species file has no abundance columns");

        _ionIndices.assign(headerCols.size(), -1);
        for (size_t i = 1; i < headerCols.size() - numTrailing; i++)
        {
            auto ion = AtomUtils::parseIon(headerCols[i]);
            _ionIndices[i] = AtomUtils::ionIndex(ion.first, ion.second);
        }
        _header = headerCols;
    }
    else if (_header != headerCols)
    {
        throw FATALERROR("Cloudy: the species file header changed between runs");
    }

    return _ionIndices;
}

////////////////////////////////////////////////////////////////////

Cloudy::Cloudy(const string& runPath, const string& inputTemplate, const CloudyConfig& config,
               CloudySpeciesHeader& species)
    : _runPath(runPath), _template(inputTemplate), _config(config), _species(species)
{}

////////////////////////////////////////////////////////////////////

void Cloudy::run(const Input& input, Output& output) const
{
    createInput(input);
    execute();
    readOutput(input, output);
}

////////////////////////////////////////////////////////////////////

void Cloudy::createInput(const Input& input) const
{
    if (!System::makeDir(_runPath)) throw FATALERROR("Cloudy: could not create the run directory " + _runPath);
    createSim(input);
    createSed(input);
}

////////////////////////////////////////////////////////////////////

void Cloudy::createSim(const Input& input) const
{
    double hden = input.hden * 1e-6;                                         // 1/m3 -> 1/cm3
    double ins = 4. * M_PI * (input.radv * _config.rad.widthv).sum() * 1e3;  // W/m2 -> erg/s/cm2

    string contents = _template;
    contents = StringUtils::replace(contents, "{hden}", StringUtils::toString(hden));
    contents = StringUtils::replace(contents, "{metal}", StringUtils::toString(input.metal));
    contents = StringUtils::replace(contents, "{ins}", StringUtils::toString(ins));

    string path = localPath("sim.in");
    std::ofstream sim = openOutput(path);
    sim.write(contents.data(), contents.size());
    closeChecked(sim, path);
}

////////////////////////////////////////////////////////////////////

void Cloudy::createSed(const Input& input) const
{
    // SKIRT has radv (4pi W/m2/m) and Cloudy wants Jlambda (~W/m2/m), so no conversion is needed.
    // Each bin is written as a flat segment between its two edges, nudged apart so that the
    // edges of neighbouring bins do not coincide.
    string path = localPath("sed.in");
    std::ofstream sed = openOutput(path);

    for (int i = 0; i < _config.rad.numBins; i++)
    {
        double left = _config.rad.edgev[i] * 0.9999;
        double right = _config.rad.edgev[i + 1] * 1.00001;
        double flux = max(input.radv[i], _config.rad.minFlux);

        sed << left << "\t" << flux << (i == 0 ? " Flambda" : "") << "\n";
        sed << right << "\t" << flux << "\n";
    }

    closeChecked(sed, path);
}

////////////////////////////////////////////////////////////////////

void Cloudy::execute() const
{
    string cmd = "cd \"" + _runPath + "\" && \"" + _config.execPath + "\" < sim.in > sim.out";
    if (System::execute(cmd) != 0) throw FATALERROR("Cloudy: the run in " + _runPath + " failed");
}

////////////////////////////////////////////////////////////////////

void Cloudy::readOutput(const Input& input, Output& output) const
{
    readTemp(output);
    readAbun(input, output);
    readOpac(output);
    readEmis(output);
    readLines(output);
}

////////////////////////////////////////////////////////////////////

void Cloudy::readTemp(Output& output) const
{
    string path = localPath("sim.ovr");
    std::ifstream ovr = openInput(path);

    string line;
    getline(ovr, line);  // header
    if (!getline(ovr, line)) throw FATALERROR("Cloudy: no data in " + path);

    output.temp = StringUtils::toDouble(splitRow(line, 2, path)[1]);
}

////////////////////////////////////////////////////////////////////

void Cloudy::readAbun(const Input& input, Output& output) const
{
    string path = localPath("sim.species");
    std::ifstream file = openInput(path);

    string header;
    string line;
    if (!getline(file, header) || !getline(file, line)) throw FATALERROR("Cloudy: no data in " + path);

    auto headerCols = StringUtils::split(header, "\t");
    auto dataCols = splitRow(line, headerCols.size(), path);
    const vector<int>& ionIndices = _species.columnsFor(headerCols);

    double hden = input.hden * 1e-6;  // 1/m3 -> 1/cm3
    for (size_t i = 0; i < ionIndices.size(); i++)
    {
        int ionIndex = ionIndices[i];
        if (ionIndex >= 0) output.abunv[ionIndex] = StringUtils::toDouble(dataCols[i]) / hden;
    }
}

////////////////////////////////////////////////////////////////////

void Cloudy::readSpectrum(const string& fileName, size_t col, Array& target) const
{
    string path = localPath(fileName);
    std::ifstream file = openInput(path);

    size_t numBins = _config.wav.numBins;
    size_t numRows = 0;
    string line;
    getline(file, line);  // header

    while (getline(file, line))
    {
        if (line.empty()) continue;
        if (numRows == numBins) throw FATALERROR("Cloudy: too many rows in " + path);

        auto cols = splitRow(line, col + 1, path);
        target[numBins - 1 - numRows] = StringUtils::toDouble(cols[col]);  // descending -> ascending
        numRows++;
    }

    if (numRows != numBins) throw FATALERROR("Cloudy: too few rows in " + path);
}

////////////////////////////////////////////////////////////////////

void Cloudy::readOpac(Output& output) const
{
    readSpectrum("sim.opac", 2, output.opacv);
    output.opacv *= 1e2;  // 1/cm -> 1/m
}

////////////////////////////////////////////////////////////////////

void Cloudy::readEmis(Output& output) const
{
    readSpectrum("sim.emis", 1, output.emisv);
    output.emisv *= 1e-1;                 // erg/s/cm3 -> W/m3
    output.emisv /= _config.wav.lambdav;  // W/m3 -> W/m3/m
}

////////////////////////////////////////////////////////////////////

void Cloudy::readLines(Output& output) const
{
    string path = localPath("sim.lines");
    std::ifstream file = openInput(path);

    string line;
    getline(file, line);  // header
    if (!getline(file, line)) throw FATALERROR("Cloudy: no data in " + path);

    // column 0 is the depth, the remaining columns are the requested lines
    auto cols = splitRow(line, _config.numLines + 1, path);

    for (int i = 0; i < _config.numLines; i++)
    {
        output.linev[i] = StringUtils::toDouble(cols[i + 1]) * 1e-7;  // ergs/s -> W
    }
}

////////////////////////////////////////////////////////////////////

string Cloudy::localPath(const string& fileName) const
{
    return StringUtils::joinPaths(_runPath, fileName);
}

////////////////////////////////////////////////////////////////////
