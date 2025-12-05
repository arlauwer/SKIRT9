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

    _space = new hnswlib::CloudySpace(_dim, Cloudy_dist);
    _hnswPath = StringUtils::joinPaths(_basePath, "hnsw.bin");
    _cloudyDir = StringUtils::joinPaths(_basePath, "cloudy");
    System::makeDir(_cloudyDir);

    // potentially load existing data
    load();
}

CloudyData CloudyWrapper::query(double hden, double metallicity, const Array& radField, double ins)
{
    if (radField.size() != cloudy::numBins) throw FATALERROR("CloudyWrapper::query: wrong number of radfield values");

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
            std::cout << "MATCH: " << label << std::endl;
        }

        return _cloudys[label];
    }
    else
    {
        _mutex.unlock();
        std::cout << "INTERPOLATE" << std::endl;
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
    // save hnsw
    if (_hnsw) _hnsw->saveIndex(StringUtils::joinPaths(_basePath, "hnsw.bin"));

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
        _hnsw = new hnswlib::HierarchicalNSW<double>(_space, _hnswPath);
    else
        _hnsw = new hnswlib::HierarchicalNSW<double>(_space, _max_elements, _M, _ef_const);

    // load cloudy
    for (auto& file : System::filesInDirectory(_cloudyDir))
    {
        int label = StringUtils::toInt(StringUtils::split(file, ".")[0]);
        int current = _next_label.load();
        while (current < label + 1 && !_next_label.compare_exchange_weak(current, label + 1))
        {
        }

        _cloudys.emplace(std::piecewise_construct, std::forward_as_tuple(label), std::forward_as_tuple());
        _dones[label] = true;
        CloudyData& data = _cloudys[label];

        string path = StringUtils::joinPaths(_cloudyDir, file);

        std::ifstream in = System::ifstream(path);
        in >> data;
        in.close();
    }
}
