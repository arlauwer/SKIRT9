#ifndef CLOUDY_WRAPPER_HPP
#define CLOUDY_WRAPPER_HPP

#include "Cloudy.hpp"
#include "NNIndex.hpp"
#include <deque>
#include <future>
#include <mutex>
#include <unordered_map>

class CloudyWrapper
{
public:
    ~CloudyWrapper();

    void setup(const CloudyConfig& config, const string& nnIndexPath);

    void save();

    void load();

    Cloudy::Output query(const Cloudy::Input& input);

private:
    void loadTemplate();

    // paths
    string _basePath;
    string _runsPath;
    string _nnIndexPath;

    // cloudy
    string _template;
    CloudyConfig _cloudyConfig;

    // cloudy data
    std::deque<Cloudy::Output> _outputs;
    Cloudy::Output _empty;

    // hnsw index
    NNIndex _nnIndex;

    // thread safe
    std::mutex _mutex;                                                // lock for query
    size_t _current_label{0};                                         // technically no atomic needed
    std::unordered_map<size_t, std::shared_future<size_t>> _pending;  // hnsw label -> future index into _outputs
};

#endif
