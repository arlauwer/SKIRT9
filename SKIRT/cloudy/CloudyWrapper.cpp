#include "CloudyWrapper.hpp"
#include "Cloudy.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "hnswlib.h"
#include "space_cloudy.h"
#include <atomic>
#include <cmath>
#include <mutex>
#include <sstream>
#include <thread>
#include <utility>

std::atomic<int> CloudyWrapper::_next_uid{0};
std::atomic<int> CloudyWrapper::_next_label{0};

namespace
{
    double hden_factor = 1e0;
    double metallicity_factor = 1e0;
    double rad_factor = 1e0;

    // cloudy distance function
    static double Cloudy_dist(const void* pVect1v, const void* pVect2v, const void* /*qty_ptr*/)
    {
        double* pVect1 = (double*)pVect1v;
        double* pVect2 = (double*)pVect2v;
        // size_t qty = *((size_t*)qty_ptr);

        double res = 0;

        // log n
        // fix this: n1 nor n2 should ever be zero!!!!
        if (*pVect1 != 0. && *pVect2 != 0.) res += hden_factor * std::abs(std::log10(*pVect1 / *pVect2));
        pVect1++;
        pVect2++;

        // linear Z
        if (*pVect1 != 0. && *pVect2 != 0.) res += metallicity_factor * std::abs(*pVect1 - *pVect2);
        pVect1++;
        pVect2++;

        // log rad
        for (size_t i = 0; i < cloudy::numBins; i++)
        {
            // should never be 0!
            res += rad_factor * std::abs(std::log10(*pVect1 / *pVect2));
            pVect1++;
            pVect2++;
        }
        return res;
    }
}

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
    if (_template.empty()) throw FATALERROR("The template file is empty");

    _space = new hnswlib::CloudySpace(_dim, Cloudy_dist);
    _hnswPath = StringUtils::joinPaths(_basePath, "hnsw.bin");

    _empty.temperature = 0.;
    _empty.abundances.resize(cloudy::numIons, 0.);
    _empty.opacities.resize(cloudy::numLambda, 0.);
    _empty.emissivities.resize(cloudy::numLambda, 0.);

    // potentially load existing data
    load();
}

CloudyData CloudyWrapper::query(double hden, double metallicity, const Array& radField, double ins)
{
    if (radField.size() != cloudy::numBins) throw FATALERROR("CloudyWrapper::query: wrong number of radfield values");

    if (hden == 0.) return _empty;

    // create query point
    vector<double> point(_dim);
    point[0] = hden;
    point[1] = metallicity;
    for (int i = 0; i < cloudy::numBins; i++) point[2 + i] = max(radField[i], cloudy::minRad);

    // there is some optimization to be done, since worst case we can have eg.
    // 8 threads, 7 are waiting to interpolate and 1 is doing perform
    // ideally we perform everything first, i.e. do a search for all cells and then determine
    // what performs have to be done

    _mutex.lock();

    // do knn search
    vector<std::pair<double, hnswlib::labeltype>> knn = _hnsw->searchKnnCloserFirst(point.data(), _k);

    // check if further than max distance
    if (knn.size() < _k || knn[_k - 1].first > _max_dist)
    {
        int label = -1;
        // check for exact match
        for (auto& nn : knn)
        {
            auto t = _hnsw->getDataByLabel<double>(nn.second);
            double dist = Cloudy_dist(point.data(), t.data(), nullptr);
            if (dist < 1e-6)  // min dist
            {
                label = nn.second;
                break;
            }
        }

        if (label == -1)
        {
            label = _next_label++;
            _hnsw->addPoint(point.data(), label);
            std::cout << "PERFORM: " << label << std::endl;

            _cloudys[label] = CloudyData();
            _dones[label] = false;

            _mutex.unlock();

            Cloudy cloudy(_next_uid++, _runsPath, _template, hden, metallicity, radField, ins);
            cloudy.perform(_cloudys[label]);
            _dones[label] = true;

            std::cout << "PERFORM DONE: " << label << std::endl;
        }
        else
        {
            _mutex.unlock();
            while (!_dones[label]) std::this_thread::sleep_for(std::chrono::milliseconds(10));
            std::cout << "MATCH: " << label << std::endl;
        }

        return _cloudys[label];
    }
    else
    {
        _mutex.unlock();
        // interpolate
        CloudyData data;

        double total_dist = 0.f;
        for (auto& pair : knn) total_dist += pair.first;

        for (auto& pair : knn)
        {
            double weight = total_dist ? (total_dist - pair.first) / total_dist : 1.;
            int id = pair.second;

            // interpolate using dist
            CloudyData& cloudy = _cloudys[id];
            while (!_dones[id]) std::this_thread::sleep_for(std::chrono::milliseconds(10));
            data.temperature += weight * cloudy.temperature;
            data.abundances += weight * cloudy.abundances;
            data.opacities += weight * cloudy.opacities;
            data.emissivities += weight * cloudy.emissivities;
        }

        return data;
    }
}

void CloudyWrapper::save()
{
    if (_hnsw) _hnsw->saveIndex(StringUtils::joinPaths(_basePath, "hnsw.bin"));

    string binPath = StringUtils::joinPaths(_basePath, "cloudy.bin");
    std::ofstream out(binPath, std::ios::binary);

    int count = _cloudys.size();
    out.write(reinterpret_cast<char*>(&count), sizeof(int));

    for (auto it = _cloudys.begin(); it != _cloudys.end(); ++it)
    {
        int label = it->first;
        CloudyData& data = it->second;

        out.write(reinterpret_cast<char*>(&label), sizeof(int));

        int nA = data.abundances.size();
        int nO = data.opacities.size();
        int nE = data.emissivities.size();

        out.write(reinterpret_cast<char*>(&nA), sizeof(int));
        out.write(reinterpret_cast<char*>(&nO), sizeof(int));
        out.write(reinterpret_cast<char*>(&nE), sizeof(int));

        out.write(reinterpret_cast<char*>(&data.temperature), sizeof(double));

        if (nA > 0) out.write(reinterpret_cast<const char*>(&data.abundances[0]), nA * sizeof(double));
        if (nO > 0) out.write(reinterpret_cast<const char*>(&data.opacities[0]), nO * sizeof(double));
        if (nE > 0) out.write(reinterpret_cast<const char*>(&data.emissivities[0]), nE * sizeof(double));
    }
}

void CloudyWrapper::load()
{
    if (System::isFile(_hnswPath))
        _hnsw = new hnswlib::HierarchicalNSW<double>(_space, _hnswPath);
    else
        _hnsw = new hnswlib::HierarchicalNSW<double>(_space, _max_elements, _M, _ef_const);

    string binPath = StringUtils::joinPaths(_basePath, "cloudy.bin");
    if (!System::isFile(binPath)) return;

    std::ifstream in(binPath, std::ios::binary);

    int count;
    in.read(reinterpret_cast<char*>(&count), sizeof(int));

    for (int i = 0; i < count; i++)
    {
        int label;
        in.read(reinterpret_cast<char*>(&label), sizeof(int));

        int nA, nO, nE;
        in.read(reinterpret_cast<char*>(&nA), sizeof(int));
        in.read(reinterpret_cast<char*>(&nO), sizeof(int));
        in.read(reinterpret_cast<char*>(&nE), sizeof(int));

        CloudyData data;
        data.abundances.resize(nA);
        data.opacities.resize(nO);
        data.emissivities.resize(nE);

        in.read(reinterpret_cast<char*>(&data.temperature), sizeof(double));

        if (nA > 0) in.read(reinterpret_cast<char*>(&data.abundances[0]), nA * sizeof(double));
        if (nO > 0) in.read(reinterpret_cast<char*>(&data.opacities[0]), nO * sizeof(double));
        if (nE > 0) in.read(reinterpret_cast<char*>(&data.emissivities[0]), nE * sizeof(double));

        _cloudys[label] = std::move(data);
        _dones[label] = true;

        int current = _next_label.load();
        while (current < label + 1 && !_next_label.compare_exchange_weak(current, label + 1))
        {
        }
    }
}
