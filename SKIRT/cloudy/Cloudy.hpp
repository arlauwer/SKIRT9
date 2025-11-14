#ifndef CLOUDY_HPP
#define CLOUDY_HPP

#include "Array.hpp"
#include "NR.hpp"
#include <array>

namespace cloudy
{
    constexpr int numBins = 4;
    constexpr std::array<double, 5> edges{1e0, 1e1, 1e2, 1e3, 1e4};
    constexpr int numIons = 465;
}

class CloudyData
{
public:
    CloudyData();

    double abundance(int Z, int N) const
    {
        int i = Z * (Z - 1) / 2 + (N - 1);
        return abundances[i];
    }

    double opacity(const Array& lambdas, double lambda) const
    {
        int l = NR::locateClip(lambdas, lambda);
        return opacities[l];
    }

    double emissivity(const Array& lambdas, double lambda) const
    {
        int l = NR::locateClip(lambdas, lambda);
        return emissivities[l];
    }

    double temperature;
    Array abundances;
    Array opacities;
    Array emissivities;
};

class Cloudy
{
public:
    Cloudy(int uid, string runsPath, const string& temp, double hden, double metallicity, const Array& radField);

    void perform();

    CloudyData& data() { return _data; }

private:
    void setup();

    void run();

    void read();

    // maybe use ions from XRayIonicGasMix and put it here for consistency?
    std::pair<int, int> getElement(string species);

    // ================== Input Data ==================
    double _hden;
    double _metallicity;
    Array _radField;

    // ================== Output Data ==================
    CloudyData _data;

    // ================== Internal Data ==================
    int _uid;
    string _path;
    const string& _template;
};

#endif
