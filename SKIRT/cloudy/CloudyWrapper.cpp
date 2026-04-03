#include "CloudyWrapper.hpp"
#include "Array.hpp"
#include "Cloudy.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <cstdint>
#include <sstream>

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

    constexpr static double fhden = 1e0;
    constexpr static double fmetal = 1e0;
    constexpr static double frad = 1e0;

    static double Cloudy_dist(const void* pVect1v, const void* pVect2v, const void* qty_ptr)
    {
        double* pVect1 = (double*)pVect1v;
        double* pVect2 = (double*)pVect2v;
        size_t qty = *((size_t*)qty_ptr);

        double res = 0;

        // log n
        // fix this: n1 nor n2 should ever be zero!!!!
        res += fhden * std::abs(std::log10(*pVect1 / *pVect2));
        pVect1++;
        pVect2++;

        // linear Z
        res += fmetal * std::abs(*pVect1 - *pVect2);
        pVect1++;
        pVect2++;

        // log rad
        for (size_t i = 0; i < qty - 2; i++)
        {
            // should never be 0!
            res += frad * std::abs(std::log10(*pVect1 / *pVect2));
            pVect1++;
            pVect2++;
        }
        return res;
    }
}

////////////////////////////////////////////////////////////////////

CloudyWrapper::~CloudyWrapper()
{
    save();
}

////////////////////////////////////////////////////////////////////

void CloudyWrapper::setup(const CloudyConfig& config, const string& nnIndexPath)
{
    _cloudyConfig = config;
    _nnIndexPath = nnIndexPath;
    _basePath = StringUtils::dirPath(nnIndexPath);
    _runsPath = StringUtils::joinPaths(_basePath, "runs");
    if (!System::makeDir(_runsPath)) throw FATALERROR("Could not create the runs directory");

    loadTemplate();

    int dim = _cloudyConfig.numDims;
    _nnIndex.setup(_nnIndexPath, dim, Cloudy_dist);

    _empty.temp = 0.;
    _empty.abunv.resize(_cloudyConfig.numIons, 0.);
    _empty.opacv.resize(_cloudyConfig.numLambda, 0.);
    _empty.emisv.resize(_cloudyConfig.numLambda + 2, 0.);

    // load existing data if available
    load();
}

////////////////////////////////////////////////////////////////////

void CloudyWrapper::save()
{
    _nnIndex.save();

    std::ofstream out(StringUtils::joinPaths(_basePath, "cloudy.out"), std::ios::binary);
    size_t count = _outputs.size();

    write(out, count);
    for (Cloudy::Output& data : _outputs)
    {
        write(out, data.temp);
        writeArray(out, data.abunv);
        writeArray(out, data.opacv);
        writeArray(out, data.emisv);
    }
}

////////////////////////////////////////////////////////////////////

void CloudyWrapper::load()
{
    _nnIndex.load();

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
    }
    _current_label = count;
}

////////////////////////////////////////////////////////////////////

Cloudy::Output CloudyWrapper::query(const Cloudy::Input& input)
{
    const auto& hden = input.hden;
    const auto& metallicity = input.metal;
    const auto& radv = input.radv;

    if (hden == 0.) return _empty;

    // create query point
    Array point(_cloudyConfig.numDims);
    point[0] = hden;
    point[1] = metallicity;
    point[std::slice(2, radv.size(), 1)] = radv;

    double* pointData = std::begin(point);

    std::unique_lock<std::mutex> lock(_mutex);  // lock

    // find approximate nearest neighbours
    vector<std::pair<double, size_t>> knn = _nnIndex.query(pointData);

    // make new cloudy point
    if (knn.empty())
    {
        size_t label = _current_label++;
        _nnIndex.addPoint(pointData, label);

        std::promise<size_t> promise;
        _pending[label] = promise.get_future().share();
        lock.unlock();  // unlock

        _outputs.emplace_back();
        Cloudy::Output& output = _outputs.back();

        string cloudyPath = "";  // TODO
        Cloudy cloudy(cloudyPath, _template, _cloudyConfig);
        cloudy.createInput(input);
        cloudy.execute();
        cloudy.readOutput(output);

        {
            std::unique_lock<std::mutex> lock(_mutex);  // lock
            promise.set_value(label);
            _pending.erase(label);
        }  // unlock

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
            label = fut.get();  // wait until Cloudy run is done
        }
        // if not pending, return
        else
        {
            lock.unlock();
        }

        return _outputs[label];
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
            if (it != _pending.end()) futs.emplace_back(it->second);
            // if not pending, add empty (invalid) future
            else
                futs.emplace_back();
        }
        lock.unlock();  // unlock

        // perform interpolation with shared futures
        Cloudy::Output interp;
        for (size_t i = 0; i < knn.size(); i++)
        {
            // get label from valid futures (wait) or from knn
            size_t label = futs[i].valid() ? futs[i].get() : knn[i].second;
            double weight = total_dist ? (total_dist - knn[i].first) / total_dist : 1.;
            const Cloudy::Output& output = _outputs[label];
            interp.temp += weight * output.temp;
            interp.abunv += weight * output.abunv;
            interp.opacv += weight * output.opacv;
            interp.emisv += weight * output.emisv;
        }
        return interp;
    }

    return _empty;
}

////////////////////////////////////////////////////////////////////

void CloudyWrapper::loadTemplate()
{
    std::ifstream in = System::ifstream(StringUtils::joinPaths(_basePath, "template.in"));
    std::ostringstream ss;
    ss << in.rdbuf();
    in.close();
    _template = ss.str();
    if (_template.empty()) throw FATALERROR("The template file is empty");
}

////////////////////////////////////////////////////////////////////
