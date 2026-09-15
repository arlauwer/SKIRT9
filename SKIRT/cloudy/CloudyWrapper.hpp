#ifndef CLOUDY_WRAPPER_HPP
#define CLOUDY_WRAPPER_HPP

#include "Cloudy.hpp"
#include <mutex>

////////////////////////////////////////////////////////////////////

class CloudyWrapper
{
public:
    void setup(CloudyConfig* config, const string& basePath);

    Cloudy::Output query(const Cloudy::Input& input);

    const Cloudy::Output& empty() const { return _empty; }

    //======== Helper Functions =======

private:
    void loadTemplate();

    const string& runPathForThisThread();

    //======================== Data Members ========================

    CloudyConfig* _config{nullptr};

    string _basePath;
    string _runsPath;
    string _template;

    Cloudy::Output _empty;
    CloudySpeciesHeader _species;

    std::mutex _mutex;  // guards _nextIndex and directory creation only
    int _nextIndex{0};
};

////////////////////////////////////////////////////////////////////

#endif
