/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Cloudy.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <fstream>

////////////////////////////////////////////////////////////////////

namespace
{
    // the files that a run writes in its directory: the input files, Cloudy's standard output, and the save files
    // requested by the input template
    const char* const runFileNames[] = {"sim.in",      "sed.in",   "sim.out",  "sim.ovr",
                                        "sim.species", "sim.emis", "sim.lines"};

    // the header line of the file written by "save diffuse continuum"; the columns are the photon energy in Ry and
    // the local continuum emission, the line emission binned onto the energy grid, and their sum, all as nu j_nu in
    // erg/s/cm3
    const char* const emisHeader = "#energy/Ryd\tConEmitLocal\tDiffuseLineEmission\tTotal";
    constexpr size_t emisTotalColumn = 3;

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
    emisv.resize(config.wav.numBins, 0.);
    linev.resize(config.numLines, 0.);
}

////////////////////////////////////////////////////////////////////

const vector<int>& CloudySpeciesHeader::columnsFor(const vector<string>& headerCols, const CloudyConfig& config)
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
            int index = config.ionIndex(headerCols[i]);
            if (index >= config.numIons) throw FATALERROR("Cloudy: ion index out of range for " + headerCols[i]);
            _ionIndices[i] = index;
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
    readOutput(output);
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
    // the intensity command sets the level of the spectrum to the total intensity 4 pi J = 4 pi sum(J_lambda dlambda);
    // it is integrated over the same floored field as written to the SED file, so that the Cloudy input depends only
    // on the floored field (as does the lookup table point), and is positive even for a zero field
    Array radv(_config.rad.numBins);
    for (int i = 0; i < _config.rad.numBins; i++) radv[i] = max(input.radv[i], _config.rad.minFlux);

    double hden = input.hden * 1e-6;                                   // 1/m3 -> 1/cm3
    double ins = 4. * M_PI * (radv * _config.rad.widthv).sum() * 1e3;  // W/m2 -> erg/s/cm2

    string contents = _template;
    contents = StringUtils::replace(contents, "{hden}", StringUtils::toString(hden));
    contents = StringUtils::replace(contents, "{metal}", StringUtils::toString(input.metal));
    contents = StringUtils::replace(contents, "{ins}", StringUtils::toString(ins));
    contents = StringUtils::replace(contents, "{minRy}", StringUtils::toString(_config.rad.borderv.min()));
    contents = StringUtils::replace(contents, "{maxRy}", StringUtils::toString(_config.rad.borderv.max()));

    string path = localPath("sim.in");
    std::ofstream sim = openOutput(path);
    sim.write(contents.data(), contents.size());
    closeChecked(sim, path);
}

////////////////////////////////////////////////////////////////////

void Cloudy::createSed(const Input& input) const
{
    // The SED file only sets the shape of the spectrum, because the intensity command in the input file sets its
    // level (see createSim()), so the mean intensity J_lambda (W/m2/m/sr) is written without conversion. Cloudy
    // expects ascending energy, so the bins are written in reverse order (they are in descending energy). Each bin
    // is written as a flat segment between its two edges, nudged inwards so that the edges of neighbouring bins do
    // not coincide.
    string path = localPath("sed.in");
    std::ofstream sed = openOutput(path);

    for (int i = _config.rad.numBins - 1; i >= 0; i--)
    {
        double low = _config.rad.borderv[i + 1] * 1.00001;
        double high = _config.rad.borderv[i] * 0.9999;
        double flux = max(input.radv[i], _config.rad.minFlux);

        sed << low << "\t" << flux << (i == _config.rad.numBins - 1 ? " Flambda" : "") << "\n";
        sed << high << "\t" << flux << "\n";
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

void Cloudy::readOutput(Output& output) const
{
    readTemp(output);
    readAbun(output);
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

void Cloudy::readAbun(Output& output) const
{
    string path = localPath("sim.species");
    std::ifstream file = openInput(path);

    string header;
    string line;
    if (!getline(file, header) || !getline(file, line)) throw FATALERROR("Cloudy: no data in " + path);

    auto headerCols = StringUtils::split(header, "\t");
    auto dataCols = splitRow(line, headerCols.size(), path);
    const vector<int>& ionIndices = _species.columnsFor(headerCols, _config);

    for (size_t i = 0; i < ionIndices.size(); i++)
    {
        int ionIndex = ionIndices[i];
        double abundance = StringUtils::toDouble(dataCols[i]) * 1e6;  // 1/cm3 -> 1/m3
        if (ionIndex >= 0) output.abunv[ionIndex] = abundance;
    }
}

////////////////////////////////////////////////////////////////////

void Cloudy::readSpectrum(const string& fileName, const string& header, size_t col, Array& target) const
{
    string path = localPath(fileName);
    std::ifstream file = openInput(path);

    size_t numBins = _config.wav.numBins;
    size_t numRows = 0;
    string line;
    getline(file, line);
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line != header) throw FATALERROR("Cloudy: unexpected header in " + path + ": " + line);

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

void Cloudy::readEmis(Output& output) const
{
    // the total column holds the continuum plus the binned lines (but no inner-shell lines)
    readSpectrum("sim.emis", emisHeader, emisTotalColumn, output.emisv);
    output.emisv *= 1e-1;                 // erg/s/cm3 -> W/m3
    output.emisv /= _config.wav.lambdav;  // W/m3 -> W/m3/m
}

////////////////////////////////////////////////////////////////////

void Cloudy::readLines(Output& output) const
{
    // without requested lines, the template has no "save lines emissivity" command
    if (_config.numLines == 0) return;

    string path = localPath("sim.lines");
    std::ifstream file = openInput(path);

    string line;
    getline(file, line);  // header
    if (!getline(file, line)) throw FATALERROR("Cloudy: no data in " + path);

    // column 0 is the depth, the remaining columns are the requested lines
    auto cols = splitRow(line, _config.numLines + 1, path);

    for (int i = 0; i < _config.numLines; i++)
    {
        output.linev[i] = StringUtils::toDouble(cols[i + 1]) * 1e-1;  // ergs/s/cm3 -> W/m3
    }
}

////////////////////////////////////////////////////////////////////

void Cloudy::removeRunFiles(const string& runPath)
{
    for (const char* name : runFileNames) System::removeFile(StringUtils::joinPaths(runPath, name));
}

////////////////////////////////////////////////////////////////////

string Cloudy::localPath(const string& fileName) const
{
    return StringUtils::joinPaths(_runPath, fileName);
}

////////////////////////////////////////////////////////////////////
