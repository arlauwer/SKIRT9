#include "CloudyWrapper.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <sstream>

////////////////////////////////////////////////////////////////////

void CloudyWrapper::setup(CloudyConfig* config, const string& basePath)
{
    _config = config;
    _basePath = basePath;
    _runsPath = StringUtils::joinPaths(_basePath, "runs");
    if (!System::makeDir(_runsPath)) throw FATALERROR("CloudyWrapper: could not create " + _runsPath);

    loadTemplate();

    _empty.resize(*_config);
}

////////////////////////////////////////////////////////////////////

Cloudy::Output CloudyWrapper::query(const Cloudy::Input& input)
{
    if (input.hden == 0.) return _empty;

    Cloudy::Output output;
    output.resize(*_config);

    Cloudy cloudy(runPathForThisThread(), _template, *_config, _species);
    cloudy.run(input, output);

    return output;
}

////////////////////////////////////////////////////////////////////

void CloudyWrapper::loadTemplate()
{
    string templatePath = FilePaths::resource("XRayCloudyGasMix_template.in");
    std::ifstream in = System::ifstream(templatePath);
    if (!in.is_open()) throw FATALERROR("CloudyWrapper: could not open the template file " + templatePath);

    std::ostringstream ss;
    ss << in.rdbuf();
    _template = ss.str();
    if (_template.empty()) throw FATALERROR("CloudyWrapper: the template file is empty " + templatePath);
}

////////////////////////////////////////////////////////////////////

const string& CloudyWrapper::runPathForThisThread()
{
    // each thread resolves and caches its own directory exactly once; after that this
    // function never touches the mutex again for that thread
    static thread_local string path;
    static thread_local bool created = false;

    if (!created)
    {
        int index;
        {
            std::unique_lock<std::mutex> lock(_mutex);
            index = _nextIndex++;
        }
        path = StringUtils::joinPaths(_runsPath, StringUtils::toString(index));
        if (!System::makeDir(path)) throw FATALERROR("CloudyWrapper: could not create the run directory " + path);
        created = true;
    }

    return path;
}

////////////////////////////////////////////////////////////////////
