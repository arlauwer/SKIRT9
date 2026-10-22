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

    <b>Direct mode</b>

    When the lookup table is disabled, every query invokes Cloudy.

    <b>Lookup table</b>

    When the lookup table is enabled, every Cloudy run is stored, and its input parameters are
    added to an approximate nearest-neighbour index (see the NNIndex class). Each input is mapped
    to a point with coordinates \f[ \left( w_n \log_{10} n_\mathrm{H},\; w_Z Z,\; \frac{w_J}{N}
    \log_{10} J_1,\; \dots,\; \frac{w_J}{N} \log_{10} J_N \right), \f] where \f$J_i\f$ is the
    radiation field in bin \f$i\f$ (clamped to the configured minimum flux) and \f$N\f$ is the
    number of radiation field bins. The distance between two inputs is the L1 distance between
    their points. With the weights used here (\f$w_n=1, w_Z=5, w_J=1\f$), a distance of one
    corresponds to a density difference of one dex, a metallicity difference of 0.2 times solar,
    or an average difference of one dex in the radiation field.

    For a query, the wrapper searches the \f$k\f$ nearest stored runs and then:
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
    different endianness. A table directory should be used by a single simulation at a time. */
class CloudyWrapper
{
public:
    /** Options controlling the lookup table. */
    struct LookupOptions
    {
        bool enabled{false};               // if false, every query invokes Cloudy
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
        must exist), and, if the lookup table is enabled, creates the table directory
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

    /** Returns the number of runs in the lookup table, or zero if the table is disabled. */
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
};

////////////////////////////////////////////////////////////////////

#endif
