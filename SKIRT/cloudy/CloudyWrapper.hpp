/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CLOUDYWRAPPER_HPP
#define CLOUDYWRAPPER_HPP

#include "Cloudy.hpp"
#include <atomic>
#include <deque>
#include <future>
#include <memory>
#include <mutex>

class NNIndex;

////////////////////////////////////////////////////////////////////

/** CloudyWrapper provides thread-safe access to Cloudy results for a given set of input
    parameters (hydrogen density, metallicity and radiation field). It manages the run directories
    and, optionally, a persistent lookup table of previous runs.

    <b>Modes</b>

    In \em direct mode, every query invokes Cloudy and nothing is stored. In \em record mode,
    every query invokes Cloudy as well, and every run is stored in the lookup table described
    below. A record run thus produces the same results as a direct run, and leaves a table that
    can later be used in lookup mode, or read by external tools. In \em lookup mode, stored runs
    are reused or interpolated where possible, and the runs that are performed are stored.

    <b>Lookup table</b>

    Every stored run is added to an approximate nearest-neighbour index (see the NNIndex class).
    Each input is mapped to a point with coordinates \f[ \left( w_n \log_{10} n_\mathrm{H},\; w_Z Z,\; \frac{w_J}{N}
    \log_{10} J_1,\; \dots,\; \frac{w_J}{N} \log_{10} J_N \right), \f] where \f$J_i\f$ is the
    radiation field in bin \f$i\f$ (clamped to the configured minimum flux) and \f$N\f$ is the
    number of radiation field bins. The distance between two inputs is the L1 distance between
    their points. With the weights used here (\f$w_n=1, w_Z=5, w_J=1\f$), a distance of one
    corresponds to a density difference of one dex, a metallicity difference of 0.2 times solar,
    or an average difference of one dex in the radiation field.

    In lookup mode, for a query, the wrapper searches the \f$k\f$ nearest stored runs and then:
    - if the nearest run lies within the match distance, returns that run;
    - otherwise, if all \f$k\f$ nearest runs lie within the interpolation distance, returns the
      inverse-distance weighted average of these runs;
    - otherwise, invokes Cloudy, stores the result, and returns it.

    To make results reusable for slightly different densities, the table stores abundances and
    opacities per unit hydrogen density, and continuum and line emissivities per unit hydrogen
    density squared. These quantities are rescaled with the density of the query. This
    approximation holds only for small differences in density, which is guaranteed by the
    distance thresholds.

    A query that would need a run already being performed by another thread waits for that run
    to complete rather than starting a duplicate run.

    <b>Multiple processes</b>

    Each process has its own wrapper, which performs its own runs. To share these runs, the client
    passes the data produced by exportNewRuns() in each process to importRuns() in all other
    processes (e.g., at the end of a medium state update), after which all wrappers hold all runs.
    The CloudyWrapper class itself does not depend on any interprocess communication mechanism.

    <b>Persistence</b>

    The table is stored in a directory. Next to the nearest-neighbour index, the directory
    contains a text file describing the configuration the table was built with (radiation field
    grid, wavelength grid, number of ions and lines, distance weights, and the full contents of the
    Cloudy input template). An existing table is only loaded if this description matches the
    current configuration exactly; otherwise setup() throws a FatalError, so that results computed
    for a different configuration are never silently reused. The lookup options (number of
    neighbours and distance thresholds) do not affect the stored data and may be changed freely
    between simulations.

    The runs are stored in a sequence of segment files. Each call to save() writes the runs added
    since the previous save to a new segment file; existing segment files are never modified. When
    a table is loaded, its segment files are memory-mapped rather than read into memory, so that
    a large table does not consume memory beyond what the operating system caches, and processes
    on the same computer share this cache. Only the runs added during the current simulation and
    the nearest-neighbour index are held in the memory of each process. The index is merely a
    cache: if it is missing or out of date, it is rebuilt from the coordinates stored with the
    runs.

    The binary files use the native byte order and are thus not portable between platforms with
    different endianness. A table directory should be used by a single simulation at a time.

    <b>Table format</b>

    A table directory contains the following files. All binary values use the native byte order
    of the machine that wrote them (little-endian on x86 and ARM computers).

    - <tt>table.txt</tt>: a text file describing the configuration, written once when the table
      is first saved. Each line ends with a newline character. Real numbers are written in
      scientific notation with 17 digits after the decimal point. The lines are, in this order:
      - <tt># Configuration of this SKIRT Cloudy lookup table; do not edit</tt>
      - <tt>version V</tt>: the format version, currently 2, equal to the version in the segment
        headers;
      - <tt>weights w_n w_Z w_J</tt>: the coordinate weights (see above);
      - <tt>ions I</tt>: the number of ion abundances in a record (495 for XRayCloudyGasMix);
      - <tt>lines L</tt>: the number of line emissivities in a record;
      - <tt>minimum-flux F</tt>: the floor for the radiation field, in W/m2/m;
      - <tt>radiation-bins N</tt>, followed by N+1 lines with the bin borders in Ry, in the order
        of the bins (for XRayCloudyGasMix, the radiation field wavelength grid in ascending
        wavelength, i.e. descending energy);
      - <tt>wavelengths W</tt>, followed by W lines with the characteristic wavelengths of the
        Cloudy output bins in m, ascending;
      - <tt>template</tt>, followed by the verbatim contents of the Cloudy input template up to
        the end of the file.
    - <tt>runs_000000.bin</tt>, <tt>runs_000001.bin</tt>, ...: the segment files, numbered
      consecutively from zero with six digits; each save adds one segment. A segment consists of
      a 64-byte header followed by its records, without padding:
      - bytes 0-15: the ASCII characters <tt>SKIRTCloudyTable</tt> (no terminating zero);
      - bytes 16-63: six unsigned 64-bit integers: the format version (2), the number of records
        R in this segment, the number of point coordinates D = N+2, the number of ions I, the
        number of wavelengths W, and the number of lines L;
      - from byte 64: R records of S = D + 1 + I + 2W + L 64-bit IEEE doubles each, so that the
        file size is 64 + 8RS bytes.
    - <tt>index.hnsw</tt>: the nearest-neighbour index in the hnswlib format, labelling each run
      with its position in the sequence of all records over all segments (in segment order),
      starting at zero. This file is a cache that is rebuilt from the records when it is missing or
      inconsistent; external tools can ignore it.
    - <tt>runs</tt> (or <tt>runs_P</tt> for process P in a multi-process run): the Cloudy run
      directories, which are not part of the table and are removed at the end of a simulation
      unless a Cloudy run failed.

    Each record holds the following values, where the offset is the index of the first value in
    the record (counted in doubles), \f$n_\mathrm{H}\f$ is the hydrogen density of the run in
    1/m3, and all quantities are in SI units:
    - offset 0: \f$w_n \log_{10} n_\mathrm{H}\f$;
    - offset 1: \f$w_Z Z\f$, with metallicity \f$Z\f$ relative to solar;
    - offsets 2 to D-1: \f$(w_J/N) \log_{10} \max(J_i, F)\f$ for the radiation bins
      \f$i=1\dots N\f$, in the order of table.txt, with \f$J_i\f$ in W/m2/m;
    - offset D: the gas temperature in K (not scaled with density);
    - offsets D+1 to D+I: the ion abundances divided by \f$n_\mathrm{H}\f$ (dimensionless), indexed
      as \f$Z(Z+1)/2+N-1\f$ for the ion with atomic number Z and N electrons (H+1, H+0, He+2, ...);
    - offsets D+1+I to D+I+W: the absorption opacities divided by \f$n_\mathrm{H}\f$, in m2, on
      the wavelengths of table.txt;
    - offsets D+1+I+W to D+I+2W: the continuum emissivities divided by \f$n_\mathrm{H}^2\f$, in
      W m3/m (i.e. (W/m3/m) / (1/m3)^2), on the same wavelengths;
    - offsets D+1+I+2W to S-1: the line emissivities divided by \f$n_\mathrm{H}^2\f$, in W m3, in
      the order of the lines requested by the template.

    The absolute quantities of a stored run are recovered by multiplying with
    \f$n_\mathrm{H} = 10^{x_0/w_n}\f$ (or its square), where \f$x_0\f$ is the record value at
    offset 0. A record mode run and a lookup mode run write identical formats. */
class CloudyWrapper
{
public:
    /** Options controlling the lookup table. */
    struct LookupOptions
    {
        /** The ways in which queries are answered: Direct invokes Cloudy for every query and stores
            nothing; Record invokes Cloudy for every query and stores every run in the table;
            Lookup reuses or interpolates stored runs where possible, as described in the class
            header, and stores the runs it performs. */
        enum class Mode { Direct, Record, Lookup };

        Mode mode{Mode::Direct};
        int numNeighbors{4};               // number of neighbours used for interpolation
        double matchDistance{0.};          // runs closer than this are reused as is
        double interpolationDistance{0.};  // all neighbours must be closer than this to interpolate
    };

    /** Counters for the different ways in which queries have been answered. */
    struct Statistics
    {
        size_t numRuns{0};
        size_t numMatches{0};
        size_t numInterpolations{0};
    };

    CloudyWrapper();

    ~CloudyWrapper();

    CloudyWrapper(const CloudyWrapper&) = delete;
    CloudyWrapper& operator=(const CloudyWrapper&) = delete;

    /** Sets up the wrapper. The function loads the Cloudy input template from \em templatePath,
        creates the directory \em runsPath for the Cloudy run directories (its parent directory
        must exist; the directory is removed again by the destructor unless a Cloudy run failed),
        and, in record or lookup mode, creates the table directory
        \em tablePath or loads the table it contains. Throws a FatalError on any failure. */
    void setup(const CloudyConfig& config, string templatePath, string tablePath, string runsPath,
               const LookupOptions& lookup);

    /** Returns the Cloudy results for the specified input, either by invoking Cloudy or from the
        lookup table as described in the class header. The radiation field must have the number of
        bins specified in the configuration. If the density is zero, the function returns empty
        results without invoking Cloudy. This function is thread-safe. */
    Cloudy::Output query(const Cloudy::Input& input);

    /** Returns empty results, i.e. results with all values set to zero. */
    const Cloudy::Output& empty() const { return _empty; }

    /** Returns true in record and lookup mode, i.e. if runs are stored in a table. */
    bool usesTable() const;

    /** Returns the number of runs in the lookup table, or zero in direct mode. */
    size_t tableSize() const;

    /** Returns the statistics accumulated since the previous call to this function, and resets
        the counters to zero. */
    Statistics takeStatistics();

    /** Appends the runs performed by this wrapper since the previous call to this function to the
        specified vector, in a form that can be passed to importRuns() of a wrapper with the same
        configuration in another process. Runs added by importRuns() are not exported again. This
        function must not be called while queries are in progress. */
    void exportNewRuns(vector<double>& data);

    /** Adds the runs exported by exportNewRuns() of a wrapper with the same configuration to the
        lookup table. This function must not be called while queries are in progress. Throws a
        FatalError if the data is inconsistent with the configuration. */
    void importRuns(const vector<double>& data);

    /** Saves the runs added to the lookup table since it was loaded or last saved, if any, to a
        new segment file in the table directory, together with the nearest-neighbour index, and
        returns true if it did so. The index is also saved if it has been rebuilt while loading. Files are written
        under a temporary name and then renamed, so that an interrupted save leaves the previous
        table intact. This function must not be called while queries are in progress. Throws a
        FatalError on any failure, including when another simulation has written to the same
        table directory in the meantime. */
    bool save();

    //======== Private helpers =======

private:
    // a run in the lookup table; the record holds the point in the nearest-neighbour index followed by the output
    // normalized to unit hydrogen density (see recordFor() for the layout), and is valid once the future is ready
    struct Entry
    {
        const double* record{nullptr};  // points into storage or into a memory-mapped segment file
        vector<double> storage;         // empty for runs loaded from a segment file
        std::promise<void> promise;
        std::shared_future<void> ready;
        bool complete{false};  // guarded by _mutex
    };

    // performs a Cloudy run in a free run directory
    void run(const Cloudy::Input& input, Cloudy::Output& output);

    // reserves and releases run directory slots
    int acquireSlot();
    void releaseSlot(int slot);

    // returns the point representing the specified input in the nearest-neighbour index
    Array pointFor(const Cloudy::Input& input) const;

    // converts between absolute and per-density quantities
    static void scale(Cloudy::Output& output, double hden);

    // stores the specified point and normalized output in a record
    vector<double> recordFor(const Array& point, const Cloudy::Output& output) const;

    // adds the weighted output stored in the specified record to the specified (resized) output
    void accumulate(Cloudy::Output& output, const double* record, double weight) const;

    // adds a complete entry with the specified record; the caller must hold _mutex
    Entry& addEntry(const double* record, vector<double>&& storage);

    // returns the text describing the configuration that determines the table contents
    string description() const;

    // loads the table from the table directory, if present
    void load();

    // writes the entries in the specified range to a segment file, or memory-maps a segment file and adds its entries
    void writeSegment(string path, size_t first, size_t last) const;
    void mapSegment(string path);

    // rebuilds the nearest-neighbour index from the points stored in the entries
    void rebuildIndex();

    string tableFile(string name) const;
    string segmentFile(int segment) const;

    //======================== Data Members ========================

private:
    CloudyConfig _config;
    LookupOptions _lookup;
    string _template;
    string _tablePath;
    string _runsPath;

    Cloudy::Output _empty;
    CloudySpeciesHeader _species;

    // run directory slots
    std::mutex _slotMutex;
    vector<char> _slotInUse;  // not vector<bool>, whose packed elements cannot be written concurrently

    // lookup table (all guarded by _mutex)
    mutable std::mutex _mutex;
    size_t _recordSize{0};  // number of values in a record
    std::unique_ptr<NNIndex> _index;
    std::deque<Entry> _entries;   // indexed by label; never invalidates references when growing
    vector<size_t> _newLabels;    // runs performed by this wrapper that have not yet been exported
    size_t _numSaved{0};          // number of entries stored in the segment files
    int _numSegments{0};          // number of segment files
    bool _indexSaved{true};       // false if the index has been rebuilt and not yet saved
    vector<string> _mappedPaths;  // memory-mapped segment files, released by the destructor

    // statistics
    std::atomic<size_t> _numRuns{0};
    std::atomic<size_t> _numMatches{0};
    std::atomic<size_t> _numInterpolations{0};

    // true if a Cloudy run failed, in which case the run directories are kept
    std::atomic<bool> _runFailed{false};
};

////////////////////////////////////////////////////////////////////

#endif
