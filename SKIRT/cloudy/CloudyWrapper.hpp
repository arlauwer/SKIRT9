#ifndef CLOUDY_WRAPPER_HPP
#define CLOUDY_WRAPPER_HPP

#include "Cloudy.hpp"
#include <atomic>
#include <mutex>
#include <unordered_map>

namespace hnswlib
{
    class CloudySpace;  // distance space for Cloudy
    template<typename T> class HierarchicalNSW;
}

class CloudyWrapper
{
public:
    ~CloudyWrapper();

    void setup(string basePath, const Array& lambda);

    CloudyData query(double hden, double metallicity, const Array& radField);

    void save();

    void load();

private:
    CloudyData perform(double hden, double metallicity, const Array& radField);

    // we need to make this thread safe somehow.
    // For a single materialmix this is easily done by making the uid atomic.
    // But for multiple mixes, not sure? Make it static?
    static std::atomic<int> _next_uid;
    static std::atomic<int> _next_label;

    string _basePath;
    string _runsPath;
    string _hnswPath;
    string _cloudyDir;

    string _template;

    // hnsw
    int _dim{2 + cloudy::numBins};
    size_t _max_elements{1000};  // set this properly!
    int _M{16};
    int _ef_const{100};
    size_t _k{1};

    hnswlib::CloudySpace* _space;
    hnswlib::HierarchicalNSW<double>* _hnsw;

    double _max_dist{1.f};

    // cloudy
    Array _lambda;
    std::unordered_map<size_t, CloudyData> _cloudys;

    std::mutex _mutex;
};

#endif
