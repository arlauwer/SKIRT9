/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "NNIndex.hpp"
#include "FatalError.hpp"
#include "System.hpp"
#include "hnswlib.h"

////////////////////////////////////////////////////////////////////

namespace
{
    // hnswlib construction parameters (see the hnswlib documentation)
    constexpr size_t numLinks = 16;         // M: number of bi-directional links per node
    constexpr size_t efConstruction = 200;  // size of the candidate list during construction
    constexpr size_t efSearch = 64;         // minimum size of the candidate list during search

    // returns the L1 distance between two points; the third argument points to the dimension
    double distanceL1(const void* a, const void* b, const void* dimPtr)
    {
        auto x = static_cast<const double*>(a);
        auto y = static_cast<const double*>(b);
        size_t numDims = *static_cast<const size_t*>(dimPtr);

        double result = 0.;
        for (size_t i = 0; i != numDims; ++i) result += std::abs(x[i] - y[i]);
        return result;
    }

    // a vector space of doubles with the L1 distance
    class SpaceL1 : public hnswlib::SpaceInterface<double>
    {
    public:
        explicit SpaceL1(size_t numDims) : _numDims(numDims) {}

        size_t get_data_size() override { return _numDims * sizeof(double); }

        hnswlib::DISTFUNC<double> get_dist_func() override { return distanceL1; }

        void* get_dist_func_param() override { return &_numDims; }

    private:
        size_t _numDims;
    };
}

////////////////////////////////////////////////////////////////////

// the hnswlib objects; the space must outlive the graph, so it is declared first
class NNIndex::Impl
{
public:
    Impl(size_t numDims, size_t capacity) : space(numDims), graph(&space, capacity, numLinks, efConstruction) {}

    Impl(size_t numDims, string path) : space(numDims), graph(&space, path) {}

    SpaceL1 space;
    hnswlib::HierarchicalNSW<double> graph;
};

////////////////////////////////////////////////////////////////////

NNIndex::NNIndex(int numDims, size_t initialCapacity)
    : _numDims(numDims), _impl(new Impl(numDims, max(initialCapacity, static_cast<size_t>(1))))
{}

////////////////////////////////////////////////////////////////////

NNIndex::~NNIndex() {}

////////////////////////////////////////////////////////////////////

size_t NNIndex::size() const
{
    return _impl->graph.getCurrentElementCount();
}

////////////////////////////////////////////////////////////////////

void NNIndex::add(const double* point, size_t label)
{
    auto& graph = _impl->graph;
    try
    {
        // grow the capacity geometrically when the index is full
        if (graph.getCurrentElementCount() >= graph.getMaxElements()) graph.resizeIndex(2 * graph.getMaxElements());
        graph.addPoint(point, label);
    }
    catch (const std::exception& error)
    {
        throw FATALERROR("Could not add a point to the nearest-neighbour index: " + string(error.what()));
    }
}

////////////////////////////////////////////////////////////////////

vector<NNIndex::Neighbor> NNIndex::search(const double* point, size_t k) const
{
    if (k == 0) return {};

    // the candidate list must be at least as long as the number of requested neighbours
    _impl->graph.setEf(max(k, efSearch));
    return _impl->graph.searchKnnCloserFirst(point, k);
}

////////////////////////////////////////////////////////////////////

void NNIndex::load(string path)
{
    if (!System::isFile(path)) throw FATALERROR("Nearest-neighbour index file does not exist: " + path);

    std::unique_ptr<Impl> impl;
    try
    {
        impl.reset(new Impl(_numDims, path));
    }
    catch (const std::exception& error)
    {
        throw FATALERROR("Could not load the nearest-neighbour index from " + path + ": " + error.what());
    }

    // hnswlib does not store the dimension, so verify the size of the stored point data instead
    const auto& graph = impl->graph;
    if (graph.label_offset_ - graph.offsetData_ != _numDims * sizeof(double))
        throw FATALERROR("The nearest-neighbour index in " + path + " has the wrong number of dimensions");

    _impl = std::move(impl);
}

////////////////////////////////////////////////////////////////////

void NNIndex::save(string path) const
{
    try
    {
        _impl->graph.saveIndex(path);
    }
    catch (const std::exception& error)
    {
        throw FATALERROR("Could not save the nearest-neighbour index to " + path + ": " + error.what());
    }
}

////////////////////////////////////////////////////////////////////
