#ifndef CLOUDY_WRAPPER_HPP
#define CLOUDY_WRAPPER_HPP

#include "Cloudy.hpp"
#include "NNIndex.hpp"
#include <deque>
#include <future>
#include <mutex>
#include <unordered_map>

// #define USE_HNSW

class CloudyWrapper
{
public:
    ~CloudyWrapper();

    void setup(const CloudyConfig& config, const string& basePath);

#ifdef USE_HNSW
    void save();

    void load();
#endif

    Cloudy::Output query(const Cloudy::Input& input);

    const Cloudy::Output& empty() const { return _empty; }

    const string& templateContent() const { return _template; }

private:
    void loadTemplate();

    int nextFreeIndex();

    // paths
    string _basePath;
    string _runsPath;

    // cloudy
    string _template;
    CloudyConfig _cloudyConfig;

    Cloudy::Output _empty;
#ifdef USE_HNSW
    // cloudy data
    std::deque<Cloudy::Output> _outputs;

    // hnsw index
    NNIndex _nnIndex;
#endif

    // thread safe
    std::mutex _mutex;                                                // lock for query
    size_t _current_label{0};                                         // technically no atomic needed
    std::unordered_map<size_t, std::shared_future<size_t>> _pending;  // hnsw label -> future index into _outputs
    vector<bool> _free;
};

#endif
