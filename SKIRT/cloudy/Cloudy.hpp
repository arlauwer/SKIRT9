#ifndef CLOUDY_HPP
#define CLOUDY_HPP

#include "Array.hpp"
#include "NR.hpp"
#include <array>
#include <iostream>

namespace cloudy
{
    constexpr int numBins = 1;
    constexpr std::array<double, 5> edges{1e0, 1e4};
    constexpr int numIons = 465;
    constexpr double minRad = 1e-10;
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

    friend std::ostream& operator<<(std::ostream& out, const CloudyData& data)
    {
        out << data.abundances.size() << '\n';
        out << data.opacities.size() << '\n';
        out << data.emissivities.size() << '\n';

        out << data.temperature << '\n';
        for (auto& abundance : data.abundances) out << abundance << '\n';
        for (auto& opacity : data.opacities) out << opacity << '\n';
        for (auto& emissivity : data.emissivities) out << emissivity << '\n';
        return out;
    }

    friend std::istream& operator>>(std::istream& in, CloudyData& data)
    {
        int numAbundances;
        int numOpacities;
        int numEmissivities;
        in >> numAbundances;
        in >> numOpacities;
        in >> numEmissivities;

        data.abundances.resize(numAbundances);
        data.opacities.resize(numOpacities);
        data.emissivities.resize(numEmissivities);

        in >> data.temperature;
        for (auto& abundance : data.abundances) in >> abundance;
        for (auto& opacity : data.opacities) in >> opacity;
        for (auto& emissivity : data.emissivities) in >> emissivity;
        return in;
    }
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
