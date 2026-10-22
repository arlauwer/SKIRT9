/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NNINDEX_HPP
#define NNINDEX_HPP

#include "Basics.hpp"
#include <memory>

////////////////////////////////////////////////////////////////////

/** NNIndex is a thin wrapper around the header-only hnswlib library, which implements approximate
    nearest-neighbour search using a hierarchical navigable small world (HNSW) graph. The wrapper
    hides the hnswlib headers from client code and adds automatic growth of the index capacity.

    The index stores points with a fixed number of dimensions, each identified by a caller-chosen
    label. Distances are measured with the L1 (Manhattan) norm, \f$d(\mathbf{x},\mathbf{y}) =
    \sum_i |x_i-y_i|\f$. A client that wishes to weigh some coordinates more heavily than others
    should scale them accordingly before handing them to the index.

    This class is <b>not</b> thread-safe. The caller must serialize all access to an instance,
    including searches. */
class NNIndex final
{
public:
    /** This type holds a search result, i.e. the distance to a point and the label of that point. */
    using Neighbor = std::pair<double, size_t>;

    /** The constructor creates an empty index for points with the specified number of dimensions.
        The initial capacity is a hint only; the index grows automatically as needed. */
    explicit NNIndex(int numDims, size_t initialCapacity = 1024);

    /** The destructor releases the hnswlib data structures. */
    ~NNIndex();

    NNIndex(const NNIndex&) = delete;
    NNIndex& operator=(const NNIndex&) = delete;

    /** This function returns the number of dimensions of the points in the index. */
    int numDims() const { return _numDims; }

    /** This function returns the number of points in the index. */
    size_t size() const;

    /** This function adds a point with the specified label to the index. The point must have
        numDims() coordinates. The label must not yet be in use. */
    void add(const double* point, size_t label);

    /** This function returns (at most) the \em k approximate nearest neighbours of the specified
        query point, sorted from closest to farthest. */
    vector<Neighbor> search(const double* point, size_t k) const;

    /** This function replaces the contents of the index by the index stored in the specified
        file, which must have been created by save() for points with the same dimension. A
        FatalError is thrown if the file cannot be read. */
    void load(string path);

    /** This function saves the index to the specified file. A FatalError is thrown if the file
        cannot be written. */
    void save(string path) const;

private:
    class Impl;

    int _numDims;
    std::unique_ptr<Impl> _impl;
};

////////////////////////////////////////////////////////////////////

#endif
