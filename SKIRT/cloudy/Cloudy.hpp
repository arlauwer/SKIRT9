#ifndef CLOUDY_HPP
#define CLOUDY_HPP

#include "Array.hpp"
#include "Basics.hpp"
#include "CloudyConfig.hpp"
#include <mutex>

////////////////////////////////////////////////////////////////////

// Stores the mapping from the Cloudy columns onto the SKIRT ion indices.
class CloudySpeciesHeader
{
public:
    /** Returns, for each column of the .species file, the SKIRT ion index, or -1 if the
        column should be ignored. */
    const vector<int>& columnsFor(const vector<string>& headerCols);

private:
    std::mutex _mutex;
    vector<string> _header;
    vector<int> _ionIndices;
};

////////////////////////////////////////////////////////////////////

/** A single Cloudy run in its own directory. Cheap to construct, not thread safe;
    create one per query, in a directory no other thread is using. */
class Cloudy
{
public:
    /** Run parameters, in SKIRT (SI) units. */
    struct Input
    {
        double hden{0.};  // 1/m3
        double metal{0.};
        Array radv;  // W/m2/m, one value per radiation bin
    };

    /** Run results, in SKIRT (SI) units. */
    struct Output
    {
        double temp{0.};  // K
        Array abunv;      // relative to hden
        Array opacv;      // 1/m,     ascending wavelength
        Array emisv;      // W/m3/m,  ascending wavelength
        Array linev;      // W

        void resize(const CloudyConfig& config);
    };

    Cloudy(const string& runPath, const string& inputTemplate, const CloudyConfig& config,
           CloudySpeciesHeader& species);

    /** Writes the input files, invokes Cloudy and parses the results into \em output,
        which must already have been resized. Throws on any failure. */
    void run(const Input& input, Output& output) const;

private:
    void createInput(const Input& input) const;
    void createSim(const Input& input) const;
    void createSed(const Input& input) const;

    void execute() const;

    void readOutput(const Input& input, Output& output) const;
    void readTemp(Output& output) const;
    void readAbun(const Input& input, Output& output) const;
    void readOpac(Output& output) const;
    void readEmis(Output& output) const;
    void readLines(Output& output) const;

    /** Reads one column of a per-wavelength Cloudy file into \em target, reversing the
        row order (Cloudy writes descending wavelength, SKIRT wants ascending). */
    void readSpectrum(const string& fileName, size_t col, Array& target) const;

    string localPath(const string& fileName) const;

    string _runPath;
    const string& _template;
    const CloudyConfig& _config;
    CloudySpeciesHeader& _species;
};

////////////////////////////////////////////////////////////////////

#endif
