#ifndef NNINDEX_HPP
#define NNINDEX_HPP

#include "Basics.hpp"

namespace hnswlib
{
    template<typename T> class SpaceInterface;
    template<typename T> class HierarchicalNSW;
}

class NNIndex final
{
public:
    using DistFunc = double (*)(const void*, const void*, const void*);

    ~NNIndex();

    void setup(const string& _hnswPath, int dim, DistFunc distFunc);

    // not thread safe!
    vector<std::pair<double, size_t>> query(const double* point);

    void addPoint(const double* point, size_t label);

    void load();

    void save();

    size_t size() const;

    size_t k() const { return _k; }

    double minDist() const { return _min_dist; }

    double maxDist() const { return _max_dist; }

private:
    string _hnswPath;
    int _dim;

    // hnswlib params
    size_t _max_elements{100000};  // set this properly!
    int _M{16};
    int _ef_const{100};
    // custom params
    size_t _k{4};
    double _min_dist{1e-4};
    double _max_dist{1.};

    // hnsw
    hnswlib::SpaceInterface<double>* _space{nullptr};
    hnswlib::HierarchicalNSW<double>* _hnsw{nullptr};
};

#endif
