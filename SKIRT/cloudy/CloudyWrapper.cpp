#include "CloudyWrapper.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "hnswlib.h"
#include <atomic>
#include <sstream>

std::atomic<int> CloudyWrapper::_next_uid{0};
std::atomic<int> CloudyWrapper::_next_label{0};

CloudyWrapper::~CloudyWrapper()
{
    save();
    delete _hnsw;
    delete _space;
}

void CloudyWrapper::setup(string basePath, const Array& lambda)
{
    _lambda = lambda;

    _basePath = basePath;
    _runsPath = StringUtils::joinPaths(_basePath, "runs");
    System::makeDir(_runsPath);

    // remove files
    auto dirs = System::dirsInDirectory(_runsPath);
    for (auto dir : dirs)
    {
        auto files = System::filesInDirectory(dir);
        for (auto file : files)
        {
            if (file == "sim.in" || file == "sim.out") System::removeFile(StringUtils::joinPaths(dir, file));
        }
    }

    std::ifstream in = System::ifstream(StringUtils::joinPaths(_basePath, "template.in"));
    std::ostringstream ss;
    ss << in.rdbuf();
    in.close();
    _template = ss.str();

    _space = new hnswlib::L2Space(_dim);
    _hnswPath = StringUtils::joinPaths(_basePath, "hnsw.bin");
    _cloudyDir = StringUtils::joinPaths(_basePath, "cloudy");
    System::makeDir(_cloudyDir);

    // potentially load existing data
    load();
}

CloudyData CloudyWrapper::query(double hden, double metallicity, const Array& radField)
{
    if (radField.size() != cloudy::numBins) throw FATALERROR("CloudyWrapper::query: wrong number of radfield values");

    // create query point
    vector<float> point(_dim);
    point[0] = hden;
    point[1] = metallicity;
    for (int i = 0; i < cloudy::numBins; i++) point[2 + i] = radField[i];

    // do knn search
    vector<std::pair<float, hnswlib::labeltype>> knn = _hnsw->searchKnnCloserFirst(point.data(), _k);

    // check if further than max distance
    if (knn.size() != _k || knn[_k - 1].first > _max_dist)
    {
        // add new point
        CloudyData data = perform(hden, metallicity, radField);
        int label = _next_label++;
        _hnsw->addPoint(point.data(), label);
        _cloudys[label] = data;

        ////// debug comparison of nn and point //////
        if (knn.size() > 0)
        {
            size_t label = knn[0].second;

            vector<float> input1 = _hnsw->getDataByLabel<float>(label);
            vector<float> input2 = point;

            auto it = _cloudys.find(label);
            if (it != _cloudys.end())
            {
                CloudyData output1 = data;
                CloudyData output2 = it->second;
            }
        }
        /////////////////////////////////////////////

        return data;
    }
    else
    {
        // nn
        size_t label = knn[0].second;
        auto it = _cloudys.find(label);
        if (it == _cloudys.end()) throw FATALERROR("CloudyWrapper::query: label not found");
        CloudyData data = it->second;

        // interpolate
        // CloudyData data;

        // float total_dist = 0.f;
        // for (auto& pair : knn) total_dist += pair.first;

        // for (auto& pair : knn)
        // {
        //     float weight = pair.first / total_dist;
        //     int id = pair.second;

        //     // interpolate using dist
        //     CloudyData& cloudy = _cloudys[id];
        //     data.temperature += weight * cloudy.temperature;
        //     data.abundances += weight * cloudy.abundances;
        //     data.opacities += weight * cloudy.opacities;
        //     data.emissivities += weight * cloudy.emissivities;
        // }
        return data;
    }
}

void CloudyWrapper::save()
{
    // save hnsw
    _hnsw->saveIndex(StringUtils::joinPaths(_basePath, "hnsw.bin"));

    // save cloudy
    for (auto& pair : _cloudys)
    {
        // save each Cloudy into its own file
        int label = pair.first;
        string path = StringUtils::joinPaths(_cloudyDir, StringUtils::toString(label));

        std::ofstream out = System::ofstream(path);
        out << pair.second;
        out.close();
    }
}

void CloudyWrapper::load()
{
    // load hnsw
    if (System::isFile(_hnswPath))
        _hnsw = new hnswlib::HierarchicalNSW<float>(_space, _hnswPath);
    else
        _hnsw = new hnswlib::HierarchicalNSW<float>(_space, _max_elements, _M, _ef_const);

    // load cloudy
    for (auto& file : System::filesInDirectory(_cloudyDir))
    {
        int label = StringUtils::toInt(StringUtils::split(file, ".")[0]);
        int current = _next_label.load();
        while (current < label + 1 && !_next_label.compare_exchange_weak(current, label + 1))
        {
        }

        CloudyData data;

        string path = StringUtils::joinPaths(_cloudyDir, file);

        std::ifstream in = System::ifstream(path);
        in >> data;
        in.close();

        _cloudys[label] = data;
    }
}

CloudyData CloudyWrapper::perform(double hden, double metallicity, const Array& radField)
{
    int uid = _next_uid++;
    Cloudy cloudy(uid, _runsPath, _template, hden, metallicity, radField);
    cloudy.perform();
    return cloudy.data();
}
