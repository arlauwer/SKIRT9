#include "NNIndex.hpp"
#include "System.hpp"
#include "hnswlib.h"

////////////////////////////////////////////////////////////////////

namespace
{
    class CustomSpace : public hnswlib::SpaceInterface<double>
    {
        size_t data_size_;
        size_t dim_;
        hnswlib::DISTFUNC<double> dist_func_;

    public:
        CustomSpace(size_t dim, NNIndex::DistFunc dist_func)
        {
            dim_ = dim;
            data_size_ = dim * sizeof(double);
            dist_func_ = dist_func;
        }

        size_t get_data_size() override { return data_size_; }

        hnswlib::DISTFUNC<double> get_dist_func() override { return dist_func_; }

        void* get_dist_func_param() override { return &dim_; }

    private:
        constexpr static double fhden = 1e0;
        constexpr static double fmetal = 1e0;
        constexpr static double frad = 1e0;
    };
}

////////////////////////////////////////////////////////////////////

NNIndex::~NNIndex()
{
    delete _space;
    delete _hnsw;
}

////////////////////////////////////////////////////////////////////

void NNIndex::setup(const string& hnswPath, int dim, DistFunc distFunc)
{
    _hnswPath = hnswPath;
    _dim = dim;
    _space = new CustomSpace(dim, distFunc);
}

////////////////////////////////////////////////////////////////////

vector<std::pair<double, size_t>> NNIndex::query(const double* point)
{
    auto knn = _hnsw->searchKnnCloserFirst(point, _k);

    // if less than _k neighbours within _max_dist
    if (knn.size() < _k || knn[_k - 1].first > _max_dist)
    {
        // found match
        if (knn.size() > 0 && knn[0].first < _min_dist)
        {
            knn.resize(1);  // only leave the matching
            return knn;
        }
        // add new point and run cloudy
        else
        {
            return {};  // else return empty list
        }
    }
    else
    {
        return knn;  // else return the k nearest neighbours
    }
}

////////////////////////////////////////////////////////////////////

void NNIndex::addPoint(const double* point, size_t label)
{
    _hnsw->addPoint(point, label);
}

////////////////////////////////////////////////////////////////////

void NNIndex::load()
{
    if (System::isFile(_hnswPath))
        _hnsw = new hnswlib::HierarchicalNSW<double>(_space, _hnswPath);
    else
        _hnsw = new hnswlib::HierarchicalNSW<double>(_space, _max_elements, _M, _ef_const);

    _hnsw->max_elements_ = _max_elements;
    _hnsw->M_ = _M;
    _hnsw->ef_construction_ = _ef_const;
}

////////////////////////////////////////////////////////////////////

void NNIndex::save()
{
    if (_hnsw) _hnsw->saveIndex(_hnswPath);
}

////////////////////////////////////////////////////////////////////

size_t NNIndex::size() const
{
    return _hnsw->getCurrentElementCount();
}

////////////////////////////////////////////////////////////////////
