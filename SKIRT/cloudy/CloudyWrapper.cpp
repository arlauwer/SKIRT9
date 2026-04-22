#include "CloudyWrapper.hpp"
#include "Array.hpp"
#include "Cloudy.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>

////////////////////////////////////////////////////////////////////

namespace
{
    template<typename T> void write(std::ofstream& out, const T& value)
    {
        out.write(reinterpret_cast<const char*>(&value), sizeof(T));
    }

    template<typename T> void read(std::ifstream& in, T& value)
    {
        in.read(reinterpret_cast<char*>(&value), sizeof(T));
    }

    void writeArray(std::ofstream& out, const Array& arr)
    {
        uint16_t n = arr.size();
        out.write(reinterpret_cast<const char*>(&n), sizeof(uint16_t));
        if (n > 0) out.write(reinterpret_cast<const char*>(&arr[0]), n * sizeof(double));
    }

    void readArray(std::ifstream& in, Array& arr)
    {
        uint16_t n;
        in.read(reinterpret_cast<char*>(&n), sizeof(uint16_t));
        arr.resize(n);
        if (n > 0) in.read(reinterpret_cast<char*>(&arr[0]), n * sizeof(double));
    }

    constexpr static double fhden = 1.;   // 1 dex for hden is too far
    constexpr static double fmetal = 5.;  // 0.2 diff in metal is too far
    constexpr static double frad = 1.;    // 1 dex for the total rad field is too far (1/N for N bins)

    static double Cloudy_dist(const void* pVect1v, const void* pVect2v, const void* qty_ptr)
    {
        double* pVect1 = (double*)pVect1v;
        double* pVect2 = (double*)pVect2v;
        size_t qty = *((size_t*)qty_ptr);

        double res = 0;

        // log n
        res += fhden * std::abs(std::log10(*pVect1 / *pVect2));
        pVect1++;
        pVect2++;

        // linear Z
        res += fmetal * std::abs(*pVect1 - *pVect2);
        pVect1++;
        pVect2++;

        // log rad
        double rad = 0.;
        for (size_t i = 0; i < qty - 2; i++)
        {
            rad += std::abs(std::log10(*pVect1 / *pVect2));
            pVect1++;
            pVect2++;
        }
        res += frad * rad / (double)(qty - 2);

        return res;
    }
}

////////////////////////////////////////////////////////////////////

CloudyWrapper::~CloudyWrapper()
{
    save();
}

////////////////////////////////////////////////////////////////////

void CloudyWrapper::setup(const CloudyConfig& config, const string& basePath)
{
    _cloudyConfig = config;
    _basePath = basePath;
    _runsPath = StringUtils::joinPaths(_basePath, "runs");
    string nnIndexPath = StringUtils::joinPaths(_basePath, "cloudy.hnsw");
    if (!System::makeDir(_runsPath)) throw FATALERROR("Could not create the runs directory");

    loadTemplate();

    int dim = _cloudyConfig.numDims;
    _nnIndex.setup(nnIndexPath, dim, Cloudy_dist);

    _empty.temp = 0.;
    _empty.abunv.resize(_cloudyConfig.numIons, 0.);
    _empty.opacv.resize(_cloudyConfig.numLambdaBins, 0.);
    _empty.emisv.resize(_cloudyConfig.numLambdaBins, 0.);
    _empty.linev.resize(_cloudyConfig.numLines, 0.);

    // load existing data if available
    load();
}

////////////////////////////////////////////////////////////////////

void CloudyWrapper::save()
{
    _nnIndex.save();

    // binary cloudy data
    std::ofstream out(StringUtils::joinPaths(_basePath, "cloudy.out"), std::ios::binary);
    size_t count = _outputs.size();

    write(out, count);
    for (Cloudy::Output& data : _outputs)
    {
        write(out, data.temp);
        writeArray(out, data.abunv);
        writeArray(out, data.opacv);
        writeArray(out, data.emisv);
        writeArray(out, data.linev);
    }

    // meta file
    string metaPath = StringUtils::joinPaths(_basePath, "cloudy.meta");
    std::ofstream meta(metaPath);
    meta << _nnIndex.size() << std::endl;
    meta << _nnIndex.k() << std::endl;
    meta << _nnIndex.minDist() << std::endl;
    meta << _nnIndex.maxDist() << std::endl;
    meta.close();
}

////////////////////////////////////////////////////////////////////

void CloudyWrapper::load()
{
    _nnIndex.load();
    _current_label = _nnIndex.size();

    // binary cloudy data
    string binPath = StringUtils::joinPaths(_basePath, "cloudy.out");
    if (!System::isFile(binPath)) return;
    std::ifstream in(binPath, std::ios::binary);

    size_t count;
    read(in, count);
    _outputs.resize(count);
    for (Cloudy::Output& data : _outputs)
    {
        read(in, data.temp);
        readArray(in, data.abunv);
        readArray(in, data.opacv);
        readArray(in, data.emisv);
        readArray(in, data.linev);
    }

    if (count != _current_label) throw FATALERROR("CloudyWrapper::load mismatch size");

    // meta file
    string metaPath = StringUtils::joinPaths(_basePath, "cloudy.meta");
    if (!System::isFile(metaPath)) return;

    std::ifstream meta(metaPath);
    size_t size, k;
    double minDist, maxDist;
    meta >> size >> k >> minDist >> maxDist;
    meta.close();

    if (size != _nnIndex.size()) throw FATALERROR("CloudyWrapper::load mismatch size");
    if (k != _nnIndex.k()) throw FATALERROR("CloudyWrapper::load mismatch k");
    if (minDist != _nnIndex.minDist()) throw FATALERROR("CloudyWrapper::load mismatch  minDist");
    if (maxDist != _nnIndex.maxDist()) throw FATALERROR("CloudyWrapper::load mismatch  maxDist");
}

////////////////////////////////////////////////////////////////////

Cloudy::Output CloudyWrapper::query(const Cloudy::Input& input)
{
    Cloudy::Output output;  // use copy to avoid race conditions
    output.resize(_cloudyConfig.numLambdaBins, _cloudyConfig.numLines);

    const auto& hden = input.hden;
    const auto& metallicity = input.metal;
    const auto& radv = input.radv;

    if (hden == 0.) return _empty;

    // create query point
    Array point(_cloudyConfig.numDims);
    point[0] = hden;
    point[1] = metallicity;
    point[std::slice(2, radv.size(), 1)] = radv;

    // clamp radiation field to minimum value
    for (int i = 0; i < (int)radv.size(); i++) point[i + 2] = std::max(radv[i], _cloudyConfig.radMin);

    double* pointData = std::begin(point);

    // ------------- LOCK -------------
    std::unique_lock<std::mutex> lock(_mutex);

    // find approximate nearest neighbours
    // the knn vector is empty if a new cloudy run is required
    vector<std::pair<double, size_t>> knn = _nnIndex.query(pointData);

    // make new cloudy point
    if (knn.empty())
    {
        // add point to hnsw
        size_t label = _current_label++;

        _nnIndex.addPoint(pointData, label);

        // make promise for result
        std::promise<size_t> promise;
        auto future = promise.get_future().share();
        _pending[label] = future;

        // index used for naming cloudy runs
        int threadIndex = nextFreeIndex();

        _outputs.emplace_back();
        // std::deque never invalidates references
        // so can be used outside the lock
        Cloudy::Output& outputRef = _outputs.back();
        outputRef.resize(_cloudyConfig.numLambdaBins, _cloudyConfig.numLines);

        lock.unlock();
        // ------------- UNLOCK -------------

        string cloudyPath = StringUtils::joinPaths(_runsPath, StringUtils::toString(threadIndex));
        Cloudy cloudy(cloudyPath, _template, _cloudyConfig);
        cloudy.createInput(input);
        cloudy.execute();
        cloudy.readOutput(outputRef);

        {
            // ------------- LOCK -------------
            std::unique_lock<std::mutex> lock(_mutex);

            // deliver promise
            promise.set_value(label);

            // remove this (shared) promise
            _pending.erase(label);

            // let other threads know the id we used is free
            _free[threadIndex] = true;

            output = outputRef;  // copy to avoid race conditions

            // ------------- UNLOCK -------------
        }

        return output;
    }
    // found exact match
    else if (knn.size() == 1)
    {
        size_t label = knn[0].second;

        auto it = _pending.find(label);
        // if pending future, wait
        if (it != _pending.end())
        {
            auto fut = it->second;
            lock.unlock();
            // ------------- UNLOCK -------------

            label = fut.get();  // wait until Cloudy run is done
        }
        // if not pending, return
        else
        {
            lock.unlock();
            // ------------- UNLOCK -------------
        }

        {
            // ------------- LOCK -------------
            std::unique_lock<std::mutex> lock(_mutex);

            output = _outputs[label];  // copy to avoid race conditions
            // ------------- UNLOCK -------------
        }

        return output;
    }
    // interpolate
    else
    {
        double total_dist = 0.;

        // Collect all pending futures if pending
        vector<std::shared_future<size_t>> futs;
        for (auto& pair : knn)
        {
            total_dist += pair.first;
            auto it = _pending.find(pair.second);
            // if pending, add future
            if (it != _pending.end())
            {
                futs.emplace_back(it->second);
            }
            // if not pending, add empty (invalid) future
            else
            {
                futs.emplace_back();
            }
        }
        lock.unlock();
        // ------------- UNLOCK -------------

        // perform interpolation with shared futures
        for (size_t i = 0; i < knn.size(); i++)
        {
            // get label from valid futures (wait) or from knn
            size_t label = futs[i].valid() ? futs[i].get() : knn[i].second;
            double weight = (total_dist - knn[i].first) / total_dist;
            const Cloudy::Output* outputRef = nullptr;
            {
                // ------------- LOCK -------------
                std::unique_lock<std::mutex> lock(_mutex);

                outputRef = &_outputs[label];  // never invalidated

                // ------------- UNLOCK -------------
            }
            output.temp += weight * outputRef->temp;
            output.abunv += weight * outputRef->abunv;
            output.opacv += weight * outputRef->opacv;
            output.emisv += weight * outputRef->emisv;
        }

        return output;
    }

    return _empty;
}

////////////////////////////////////////////////////////////////////

void CloudyWrapper::loadTemplate()
{
    std::ifstream in = System::ifstream(StringUtils::joinPaths(_basePath, "XRayIonicGasMix_template.in"));
    std::ostringstream ss;
    ss << in.rdbuf();
    in.close();
    _template = ss.str();
    if (_template.empty()) throw FATALERROR("The template file is empty");
}

////////////////////////////////////////////////////////////////////

int CloudyWrapper::nextFreeIndex()
{
    // find first free index
    auto it = std::find(_free.begin(), _free.end(), true);
    int index = it - _free.begin();

    // if no free index found, add one and use it
    if (it == _free.end())
    {
        _free.resize(_free.size() + 1, false);
    }
    // else use the found index
    else
    {
        _free[index] = false;
    }

    return index;
}

////////////////////////////////////////////////////////////////////
