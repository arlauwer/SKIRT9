#include "CloudyWrapper.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "hnswlib.h"
#include <sstream>

std::atomic<int> CloudyWrapper::_next_uid{0};

CloudyWrapper::~CloudyWrapper()
{
    delete _hnsw;
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

    _points = new float*[_max_elements];
    _space = new hnswlib::L2Space(_dim);
    _hnsw = new hnswlib::HierarchicalNSW<float>(_space, _max_elements, _M, _ef_const);
}

CloudyData CloudyWrapper::query(double hden, double metallicity, const Array& radField)
{
    _points[_next_uid] = new float[_dim];
    _points[_next_uid][0] = hden;
    _points[_next_uid][1] = metallicity;

    if (radField.size() != cloudy::numBins) throw FATALERROR("CloudyWrapper::query: wrong number of radfield values");
    for (int i = 0; i < cloudy::numBins; i++) _points[_next_uid][2 + i] = radField[i];

    vector<std::pair<float, hnswlib::labeltype>> knn =
        _hnsw->searchKnnCloserFirst(_points[_next_uid], _k);  // problem here: segfault at 2nd iteration!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // check if further than max distance
    if (knn.size() != _k || knn[_k - 1].first > _max_dist)
    {
        Cloudy cloudy = perform(hden, metallicity, radField);
        _hnsw->addPoint(_points[_next_uid], _hnsw->data_size_);
        _cloudys.push_back(cloudy);
        return cloudy.data();
    }
    else
    {
        CloudyData data;

        float total_dist = 0.f;
        for (auto& pair : knn) total_dist += pair.first;

        for (auto& pair : knn)
        {
            float weight = pair.first / total_dist;
            int id = pair.second;

            // interpolate using dist
            CloudyData& cloudy = _cloudys[id].data();
            data.temperature += weight * cloudy.temperature;
            data.abundances += weight * cloudy.abundances;
            data.opacities += weight * cloudy.opacities;
            data.emissivities += weight * cloudy.emissivities;
        }
        return data;
    }
}

Cloudy CloudyWrapper::perform(double hden, double metallicity, const Array& radField)
{
    int uid = _next_uid++;
    Cloudy cloudy(uid, _runsPath, _template, hden, metallicity, radField);
    cloudy.perform();
    return cloudy;
}
