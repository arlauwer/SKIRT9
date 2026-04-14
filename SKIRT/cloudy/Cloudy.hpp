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
    int numLambdaBins{-1};
    Array lambdaBorderv;  // meter in ascending order

    // --- Lines ---
    int numLines{-1};
    Array lineEmisCenterv;  // m
    Array lineMassv;        // amu

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
        void resize(int numBins, int numLines)
        {
            abunv.resize(CloudyConfig::numIons, 0.);
            opacv.resize(numBins, 0.);
            emisv.resize(numBins, 0.);
            linev.resize(numLines, 0.);
        }

        double temp;
        Array abunv;  // (1/m3)
        Array opacv;  // (1/m)  ascending wavelength
        Array emisv;  // (W/m3) ascending wavelength
        Array linev;  // (W/m3)
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

    void binSegments(const string& fileName, Array& binnedSpectrum, size_t col) const;

    void readOpac(Output& output) const;

    void readEmis(Output& output) const;

    void readLines(Output& output) const;

private:
    string localPath(const string& filename) const;

    string _basePath;
    const string& _inputTemplate;
    const CloudyConfig& _cloudyConfig;
};

#endif
