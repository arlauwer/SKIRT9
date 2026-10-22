/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CLOUDY_HPP
#define CLOUDY_HPP

#include "Array.hpp"
#include "Basics.hpp"
#include "CloudyConfig.hpp"
#include <mutex>

////////////////////////////////////////////////////////////////////

/** CloudySpeciesHeader maps the columns of the Cloudy species file onto ion indices. The mapping
    is resolved once, on the first run, and shared by all subsequent runs. This class is
    thread-safe. */
class CloudySpeciesHeader
{
public:
    /** Returns, for each column of the species file with the specified header, the ion index, or
        -1 if the column should be ignored. Throws a FatalError if the header differs from the one
        seen on a previous call. */
    const vector<int>& columnsFor(const vector<string>& headerCols, const CloudyConfig& config);

private:
    std::mutex _mutex;
    vector<string> _header;
    vector<int> _ionIndices;
};

////////////////////////////////////////////////////////////////////

/** Cloudy represents a single one-zone Cloudy run in its own directory. It writes the input files,
    invokes the Cloudy executable, and parses the results, converting everything to SKIRT (SI)
    units. A Cloudy instance is cheap to construct but not thread-safe; create one per run, in a
    directory that no other thread is using at the same time. */
class Cloudy
{
public:
    /** Run parameters, in SKIRT (SI) units. */
    struct Input
    {
        double hden{0.};   // 1/m3
        double metal{0.};  // relative to solar
        Array radv;        // W/m2/m, one value per radiation bin
    };

    /** Run results, in SKIRT (SI) units. */
    struct Output
    {
        double temp{0.};  // K
        Array abunv;      // 1/m3, one value per ion
        Array opacv;      // 1/m, ascending wavelength
        Array emisv;      // W/m3/m, ascending wavelength
        Array linev;      // W/m3, one value per line

        /** Resizes all arrays according to the specified configuration and sets all values to
            zero. */
        void resize(const CloudyConfig& config);
    };

    /** Constructs a run in the directory \em runPath with the specified template contents,
        configuration and species mapping. The referenced objects must outlive this instance. */
    Cloudy(const string& runPath, const string& inputTemplate, const CloudyConfig& config,
           CloudySpeciesHeader& species);

    /** Writes the input files, invokes Cloudy and parses the results into \em output, which must
        already have been resized. Throws a FatalError on any failure. */
    void run(const Input& input, Output& output) const;

    /** Removes the files that a run writes in the directory \em runPath (input, output and save
        files, as named in the input template), and no other files. Missing files are ignored. */
    static void removeRunFiles(const string& runPath);

private:
    void createInput(const Input& input) const;
    void createSim(const Input& input) const;
    void createSed(const Input& input) const;

    void execute() const;

    void readOutput(Output& output) const;
    void readTemp(Output& output) const;
    void readAbun(Output& output) const;
    void readOpac(Output& output) const;
    void readEmis(Output& output) const;
    void readLines(Output& output) const;

    /** Reads one column of a per-wavelength Cloudy file into \em target, reversing the row order
        (Cloudy writes descending wavelength, SKIRT wants ascending). */
    void readSpectrum(const string& fileName, size_t col, Array& target) const;

    string localPath(const string& fileName) const;

    string _runPath;
    const string& _template;
    const CloudyConfig& _config;
    CloudySpeciesHeader& _species;
};

////////////////////////////////////////////////////////////////////

#endif
