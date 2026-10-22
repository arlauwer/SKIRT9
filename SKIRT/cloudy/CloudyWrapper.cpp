/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CloudyWrapper.hpp"
#include "FatalError.hpp"
#include "NNIndex.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>

////////////////////////////////////////////////////////////////////

namespace
{
    // weights of the coordinates in the nearest-neighbour index (see the class header)
    constexpr double hdenWeight = 1.;   // per dex in hydrogen density
    constexpr double metalWeight = 5.;  // per unit of (solar) metallicity
    constexpr double radWeight = 1.;    // per dex in the radiation field, averaged over all bins

    // file names in the table directory
    const char* const descriptionFileName = "table.txt";
    const char* const indexFileName = "index.hnsw";
    const char* const segmentFilePrefix = "runs_";
    const char* const segmentFileSuffix = ".bin";

    // header of a segment file: a 16-character identifier followed by 64-bit integers, so that the records that follow
    // are aligned on a multiple of 8 bytes
    const string segmentMagic = "SKIRTCloudyTable";
    constexpr uint64_t segmentVersion = 2;
    constexpr size_t numHeaderValues = 6;  // version, number of records, dims, ions, wavelengths, lines
    constexpr size_t headerSize = 16 + numHeaderValues * sizeof(uint64_t);

    // returns the number of coordinates of a point in the nearest-neighbour index
    int numDimsFor(const CloudyConfig& config)
    {
        return 2 + config.rad.numBins;
    }

    // returns the full-precision representation of a double
    string toExactString(double value)
    {
        return StringUtils::toString(value, 'e', 17);
    }

    // returns the full contents of a text file, or throws
    string readTextFile(const string& path)
    {
        std::ifstream in = System::ifstream(path);
        if (!in.is_open()) throw FATALERROR("Could not open " + path);
        std::ostringstream contents;
        contents << in.rdbuf();
        return contents.str();
    }

    // writes the full contents of a text file, or throws
    void writeTextFile(const string& path, const string& contents)
    {
        std::ofstream out = System::ofstream(path);
        out << contents;
        out.close();
        if (out.fail()) throw FATALERROR("Could not write " + path);
    }

    // renames a file, replacing any existing file with the new name, or throws
    void renameFile(const string& from, const string& to)
    {
        if (std::rename(from.c_str(), to.c_str()) == 0) return;

        // on Windows, rename does not replace an existing file
        std::remove(to.c_str());
        if (std::rename(from.c_str(), to.c_str()) != 0) throw FATALERROR("Could not rename " + from + " to " + to);
    }

    template<typename T> void writeValue(std::ofstream& out, T value)
    {
        out.write(reinterpret_cast<const char*>(&value), sizeof(T));
    }

    void writeValues(std::ofstream& out, const double* values, size_t count)
    {
        if (count) out.write(reinterpret_cast<const char*>(values), count * sizeof(double));
    }
}

////////////////////////////////////////////////////////////////////

CloudyWrapper::CloudyWrapper() {}

////////////////////////////////////////////////////////////////////

CloudyWrapper::~CloudyWrapper()
{
    for (const string& path : _mappedPaths) System::releaseMemoryMap(path);
}

////////////////////////////////////////////////////////////////////

void CloudyWrapper::setup(const CloudyConfig& config, string templatePath, string tablePath, string runsPath,
                          const LookupOptions& lookup)
{
    _config = config;
    _lookup = lookup;
    _tablePath = tablePath;
    _runsPath = runsPath;

    if (!_config.ionIndex) throw FATALERROR("CloudyWrapper: no ion mapping has been configured");
    if (_lookup.enabled)
    {
        if (_lookup.numNeighbors < 1) throw FATALERROR("CloudyWrapper: the number of neighbours must be positive");
        if (_lookup.matchDistance < 0. || _lookup.interpolationDistance < _lookup.matchDistance)
            throw FATALERROR("CloudyWrapper: the distances must satisfy 0 <= match <= interpolation");
    }

    // load the input template
    _template = readTextFile(templatePath);
    if (_template.empty()) throw FATALERROR("CloudyWrapper: the template file is empty: " + templatePath);

    // create the directory holding the run directories
    if (!System::makeDir(_runsPath)) throw FATALERROR("CloudyWrapper: could not create " + _runsPath);

    _empty.resize(_config);

    // create or load the lookup table
    if (_lookup.enabled)
    {
        if (!System::makeDir(_tablePath)) throw FATALERROR("CloudyWrapper: could not create " + _tablePath);
        _recordSize = numDimsFor(_config) + 1 + _config.numIons + 2 * _config.wav.numBins + _config.numLines;
        _index.reset(new NNIndex(numDimsFor(_config)));
        load();
    }
}

////////////////////////////////////////////////////////////////////

Cloudy::Output CloudyWrapper::query(const Cloudy::Input& input)
{
    if (input.hden <= 0.) return _empty;
    if (static_cast<int>(input.radv.size()) != _config.rad.numBins)
        throw FATALERROR("CloudyWrapper: the radiation field has the wrong number of bins");

    // direct mode: always run Cloudy
    if (!_lookup.enabled)
    {
        Cloudy::Output output;
        run(input, output);
        return output;
    }

    Array point = pointFor(input);
    std::unique_lock<std::mutex> lock(_mutex);
    auto neighbors = _index->search(&point[0], _lookup.numNeighbors);

    // reuse the nearest run if it is close enough; a distance so small that its inverse overflows is treated as a
    // match as well, so that the interpolation weights below are always finite
    if (!neighbors.empty() && (neighbors[0].first <= _lookup.matchDistance || !std::isfinite(1. / neighbors[0].first)))
    {
        const Entry& entry = _entries[neighbors[0].second];
        auto ready = entry.ready;
        lock.unlock();

        ready.get();  // wait for the run to complete if it is still in progress; rethrows on failure
        Cloudy::Output output;
        output.resize(_config);
        accumulate(output, entry.record, 1.);
        scale(output, input.hden);
        _numMatches++;
        return output;
    }

    // interpolate if all neighbours are close enough
    if (static_cast<int>(neighbors.size()) == _lookup.numNeighbors
        && neighbors.back().first <= _lookup.interpolationDistance)
    {
        vector<const Entry*> entries;
        vector<std::shared_future<void>> readies;
        for (const auto& neighbor : neighbors)
        {
            entries.push_back(&_entries[neighbor.second]);
            readies.push_back(entries.back()->ready);
        }
        lock.unlock();

        // inverse-distance weighting; all distances exceed the match distance and have a finite inverse (see above)
        Cloudy::Output output;
        output.resize(_config);
        double sumWeights = 0.;
        for (size_t i = 0; i != neighbors.size(); ++i)
        {
            readies[i].get();
            double weight = 1. / neighbors[i].first;
            accumulate(output, entries[i]->record, weight);
            sumWeights += weight;
        }
        if (!(sumWeights > 0.) || !std::isfinite(sumWeights))
            throw FATALERROR("CloudyWrapper: degenerate interpolation weights");

        output.temp /= sumWeights;
        output.abunv /= sumWeights;
        output.opacv /= sumWeights;
        output.emisv /= sumWeights;
        output.linev /= sumWeights;
        scale(output, input.hden);
        _numInterpolations++;
        return output;
    }

    // otherwise add a new entry and run Cloudy; other threads querying nearby points wait for it
    size_t label = _entries.size();
    _entries.emplace_back();
    Entry& entry = _entries.back();
    entry.ready = entry.promise.get_future().share();
    _index->add(&point[0], label);
    lock.unlock();

    // the promise must be fulfilled on every path, because threads waiting for this run would otherwise never wake up
    Cloudy::Output output;
    try
    {
        run(input, output);
        Cloudy::Output normalized = output;
        scale(normalized, 1. / input.hden);
        entry.storage = recordFor(point, normalized);
        entry.record = entry.storage.data();
    }
    catch (...)
    {
        entry.promise.set_exception(std::current_exception());
        throw;
    }
    entry.promise.set_value();

    lock.lock();
    entry.complete = true;
    _newLabels.push_back(label);
    return output;
}

////////////////////////////////////////////////////////////////////

size_t CloudyWrapper::tableSize() const
{
    std::unique_lock<std::mutex> lock(_mutex);
    return _entries.size();
}

////////////////////////////////////////////////////////////////////

CloudyWrapper::Statistics CloudyWrapper::takeStatistics()
{
    Statistics statistics;
    statistics.numRuns = _numRuns.exchange(0);
    statistics.numMatches = _numMatches.exchange(0);
    statistics.numInterpolations = _numInterpolations.exchange(0);
    return statistics;
}

////////////////////////////////////////////////////////////////////

void CloudyWrapper::exportNewRuns(vector<double>& data)
{
    std::unique_lock<std::mutex> lock(_mutex);
    if (!_lookup.enabled) return;

    data.reserve(data.size() + _newLabels.size() * _recordSize);
    for (size_t label : _newLabels)
    {
        const Entry& entry = _entries[label];
        if (!entry.complete) throw FATALERROR("CloudyWrapper: cannot export runs while Cloudy runs are in progress");
        data.insert(data.end(), entry.record, entry.record + _recordSize);
    }
    _newLabels.clear();
}

////////////////////////////////////////////////////////////////////

void CloudyWrapper::importRuns(const vector<double>& data)
{
    std::unique_lock<std::mutex> lock(_mutex);
    if (!_lookup.enabled) return;
    if (data.size() % _recordSize) throw FATALERROR("CloudyWrapper: imported runs have an unexpected size");

    for (size_t first = 0; first != data.size(); first += _recordSize)
    {
        vector<double> storage(data.begin() + first, data.begin() + first + _recordSize);
        Entry& entry = addEntry(storage.data(), std::move(storage));
        _index->add(entry.record, _entries.size() - 1);
    }
}

////////////////////////////////////////////////////////////////////

bool CloudyWrapper::save()
{
    std::unique_lock<std::mutex> lock(_mutex);
    if (!_lookup.enabled || (_numSaved == _entries.size() && _indexSaved)) return false;

    for (size_t label = _numSaved; label != _entries.size(); ++label)
        if (!_entries[label].complete)
            throw FATALERROR("CloudyWrapper: cannot save the table while Cloudy runs are in progress");

    // the description never changes for a given table, so write it only once, before any segment file
    string suffix = ".tmp";
    string descriptionPath = tableFile(descriptionFileName);
    if (!System::isFile(descriptionPath))
    {
        writeTextFile(descriptionPath + suffix, description());
        renameFile(descriptionPath + suffix, descriptionPath);
    }

    // write the new runs to a new segment file; refuse to overwrite a segment written by another simulation
    if (_numSaved != _entries.size())
    {
        string segmentPath = segmentFile(_numSegments);
        if (System::isFile(segmentPath))
            throw FATALERROR("CloudyWrapper: the lookup table in " + _tablePath
                             + " has been modified by another simulation; use a separate table directory");
        writeSegment(segmentPath + suffix, _numSaved, _entries.size());
        renameFile(segmentPath + suffix, segmentPath);
        _numSegments++;
        _numSaved = _entries.size();
    }

    // replace the index, which is rebuilt when loading the table if this step is interrupted
    string indexPath = tableFile(indexFileName);
    _index->save(indexPath + suffix);
    renameFile(indexPath + suffix, indexPath);
    _indexSaved = true;
    return true;
}

////////////////////////////////////////////////////////////////////

void CloudyWrapper::run(const Cloudy::Input& input, Cloudy::Output& output)
{
    output.resize(_config);

    // release the slot when leaving this function, also when Cloudy throws
    struct SlotGuard
    {
        CloudyWrapper* wrapper;
        int slot;
        ~SlotGuard() { wrapper->releaseSlot(slot); }
    } guard{this, acquireSlot()};

    string runPath = StringUtils::joinPaths(_runsPath, StringUtils::toString(guard.slot));
    Cloudy cloudy(runPath, _template, _config, _species);
    cloudy.run(input, output);
    _numRuns++;
}

////////////////////////////////////////////////////////////////////

int CloudyWrapper::acquireSlot()
{
    std::unique_lock<std::mutex> lock(_slotMutex);
    for (size_t slot = 0; slot != _slotInUse.size(); ++slot)
    {
        if (!_slotInUse[slot])
        {
            _slotInUse[slot] = true;
            return slot;
        }
    }
    _slotInUse.push_back(true);
    return _slotInUse.size() - 1;
}

////////////////////////////////////////////////////////////////////

void CloudyWrapper::releaseSlot(int slot)
{
    std::unique_lock<std::mutex> lock(_slotMutex);
    _slotInUse[slot] = false;
}

////////////////////////////////////////////////////////////////////

Array CloudyWrapper::pointFor(const Cloudy::Input& input) const
{
    int numBins = _config.rad.numBins;
    Array point(numDimsFor(_config));
    point[0] = hdenWeight * std::log10(input.hden);
    point[1] = metalWeight * input.metal;
    for (int i = 0; i != numBins; ++i)
        point[2 + i] = radWeight / numBins * std::log10(max(input.radv[i], _config.rad.minFlux));
    return point;
}

////////////////////////////////////////////////////////////////////

void CloudyWrapper::scale(Cloudy::Output& output, double hden)
{
    output.abunv *= hden;
    output.opacv *= hden;
    output.emisv *= hden * hden;
    output.linev *= hden * hden;
}

////////////////////////////////////////////////////////////////////

vector<double> CloudyWrapper::recordFor(const Array& point, const Cloudy::Output& output) const
{
    // layout: point, temperature, abundances, opacities, emissivities, line emissivities
    vector<double> record;
    record.reserve(_recordSize);
    record.insert(record.end(), std::begin(point), std::end(point));
    record.push_back(output.temp);
    record.insert(record.end(), std::begin(output.abunv), std::end(output.abunv));
    record.insert(record.end(), std::begin(output.opacv), std::end(output.opacv));
    record.insert(record.end(), std::begin(output.emisv), std::end(output.emisv));
    record.insert(record.end(), std::begin(output.linev), std::end(output.linev));
    return record;
}

////////////////////////////////////////////////////////////////////

void CloudyWrapper::accumulate(Cloudy::Output& output, const double* record, double weight) const
{
    const double* value = record + numDimsFor(_config);
    output.temp += weight * *value++;
    for (double& v : output.abunv) v += weight * *value++;
    for (double& v : output.opacv) v += weight * *value++;
    for (double& v : output.emisv) v += weight * *value++;
    for (double& v : output.linev) v += weight * *value++;
}

////////////////////////////////////////////////////////////////////

CloudyWrapper::Entry& CloudyWrapper::addEntry(const double* record, vector<double>&& storage)
{
    _entries.emplace_back();
    Entry& entry = _entries.back();
    entry.storage = std::move(storage);
    entry.record = record;
    entry.ready = entry.promise.get_future().share();
    entry.promise.set_value();
    entry.complete = true;
    return entry;
}

////////////////////////////////////////////////////////////////////

string CloudyWrapper::description() const
{
    std::ostringstream out;
    out << "# Configuration of this SKIRT Cloudy lookup table; do not edit\n";
    out << "version " << segmentVersion << "\n";
    out << "weights " << toExactString(hdenWeight) << " " << toExactString(metalWeight) << " "
        << toExactString(radWeight) << "\n";
    out << "ions " << _config.numIons << "\n";
    out << "lines " << _config.numLines << "\n";
    out << "minimum-flux " << toExactString(_config.rad.minFlux) << "\n";
    out << "radiation-bins " << _config.rad.numBins << "\n";
    for (double border : _config.rad.borderv) out << toExactString(border) << "\n";
    out << "wavelengths " << _config.wav.numBins << "\n";
    for (double lambda : _config.wav.lambdav) out << toExactString(lambda) << "\n";
    out << "template\n" << _template;
    return out.str();
}

////////////////////////////////////////////////////////////////////

void CloudyWrapper::load()
{
    string descriptionPath = tableFile(descriptionFileName);
    string indexPath = tableFile(indexFileName);

    if (!System::isFile(descriptionPath))
    {
        if (System::isFile(indexPath) || System::isFile(segmentFile(0)))
            throw FATALERROR("CloudyWrapper: incomplete lookup table in " + _tablePath);
        return;  // start with an empty table
    }

    if (readTextFile(descriptionPath) != description())
        throw FATALERROR("CloudyWrapper: the lookup table in " + _tablePath
                         + " was built with a different configuration; specify another table directory");

    while (System::isFile(segmentFile(_numSegments))) mapSegment(segmentFile(_numSegments++));
    _numSaved = _entries.size();

    // use the saved index if it is consistent with the runs, and otherwise rebuild it
    bool indexLoaded = false;
    if (System::isFile(indexPath))
    {
        try
        {
            _index->load(indexPath);
            indexLoaded = _index->size() == _entries.size();
        }
        catch (FatalError&)
        {
            // a damaged index is rebuilt below
        }
    }
    if (!indexLoaded) rebuildIndex();
}

////////////////////////////////////////////////////////////////////

void CloudyWrapper::writeSegment(string path, size_t first, size_t last) const
{
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out.is_open()) throw FATALERROR("CloudyWrapper: could not open " + path);

    out.write(segmentMagic.data(), segmentMagic.size());
    writeValue<uint64_t>(out, segmentVersion);
    writeValue<uint64_t>(out, last - first);
    writeValue<uint64_t>(out, numDimsFor(_config));
    writeValue<uint64_t>(out, _config.numIons);
    writeValue<uint64_t>(out, _config.wav.numBins);
    writeValue<uint64_t>(out, _config.numLines);

    for (size_t label = first; label != last; ++label) writeValues(out, _entries[label].record, _recordSize);

    out.close();
    if (out.fail()) throw FATALERROR("CloudyWrapper: could not write " + path);
}

////////////////////////////////////////////////////////////////////

void CloudyWrapper::mapSegment(string path)
{
    auto map = System::acquireMemoryMap(path);
    if (!map.first) throw FATALERROR("CloudyWrapper: could not memory-map " + path);
    _mappedPaths.push_back(path);

    const char* bytes = static_cast<const char*>(map.first);
    if (map.second < headerSize || string(bytes, segmentMagic.size()) != segmentMagic)
        throw FATALERROR("CloudyWrapper: unexpected header in " + path);

    uint64_t header[numHeaderValues];
    std::memcpy(header, bytes + segmentMagic.size(), sizeof(header));
    uint64_t count = header[1];
    if (header[0] != segmentVersion || header[2] != static_cast<uint64_t>(numDimsFor(_config))
        || header[3] != static_cast<uint64_t>(_config.numIons)
        || header[4] != static_cast<uint64_t>(_config.wav.numBins)
        || header[5] != static_cast<uint64_t>(_config.numLines))
        throw FATALERROR("CloudyWrapper: unexpected header in " + path);
    if (map.second != headerSize + count * _recordSize * sizeof(double))
        throw FATALERROR("CloudyWrapper: unexpected size of " + path);

    // the map is page-aligned and the header size is a multiple of 8 bytes, so the records are properly aligned
    const double* records = reinterpret_cast<const double*>(bytes + headerSize);
    for (uint64_t i = 0; i != count; ++i) addEntry(records + i * _recordSize, vector<double>());
}

////////////////////////////////////////////////////////////////////

void CloudyWrapper::rebuildIndex()
{
    _index.reset(new NNIndex(numDimsFor(_config), max(_entries.size(), static_cast<size_t>(1024))));
    for (size_t label = 0; label != _entries.size(); ++label) _index->add(_entries[label].record, label);
    _indexSaved = false;
}

////////////////////////////////////////////////////////////////////

string CloudyWrapper::tableFile(string name) const
{
    return StringUtils::joinPaths(_tablePath, name);
}

////////////////////////////////////////////////////////////////////

string CloudyWrapper::segmentFile(int segment) const
{
    string number = StringUtils::toString(segment);
    if (number.size() < 6) number.insert(0, 6 - number.size(), '0');
    return tableFile(segmentFilePrefix + number + segmentFileSuffix);
}

////////////////////////////////////////////////////////////////////
