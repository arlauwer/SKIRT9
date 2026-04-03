#ifndef CLOUDY_HPP
#define CLOUDY_HPP

#include "Array.hpp"
#include "Basics.hpp"

struct CloudyConfig
{
    // --- Radiation field ---
    size_t numBins;
    Array radEdges;  // rydberg
    Array radWidth;  // meter
    double radMin;   // W/m2/m

    // --- Optical properties ---
    Array lambdav;  // meter
    int numLambda;

    // --- Ions ---
    constexpr static int numIons = 465;  // might have variable number of ions later
    constexpr static int numAtoms = 30;

    // --- Cloudy ---
    string cloudyExecPath;
    size_t numDims;  // total number of parameters
};

class Cloudy
{
public:
    // using SKIRT (SI) units
    struct Input
    {
        double hden{0};
        double metal{0};
        Array radv;
    };

    // using SKIRT (SI) units
    struct Output
    {
        double temp{0};
        Array abunv;
        Array opacv;
        Array emisv;
    };

    Cloudy(string basePath, const string& inputTemplate, const CloudyConfig& config);

    void createInput(const Input& input) const;

    bool execute();

    void readOutput(Output& output) const;

private:
    void createSim(const Input& input) const;

    void createSed(const Input& input) const;

    void readTemp(Output& output) const;

    void readAbun(Output& output) const;

    void readOpac(Output& output) const;

    void readEmis(Output& output) const;

private:
    string localPath(const string& filename) const;

    string _basePath;
    const string& _inputTemplate;
    const CloudyConfig& _cloudyConfig;
};

#endif
