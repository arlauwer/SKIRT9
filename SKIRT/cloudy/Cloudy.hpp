#ifndef CLOUDY_HPP
#define CLOUDY_HPP

#include "Array.hpp"
#include "Constants.hpp"
#include "NR.hpp"
#include <array>
#include <atomic>
#include <iostream>

namespace cloudy
{
    constexpr int numBins = 1;
    constexpr std::array<double, numBins + 1> edges{1e4, 1e0};
    constexpr int numIons = 465;
    constexpr double minRad = 1e-10;
    constexpr double stc = 1e3;                // W/m2 -> erg/s/cm2
    constexpr double Rtm = Constants::iRyd();  // Ryd -> m
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
    std::atomic<bool> done;

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
    Cloudy(int uid, string runsPath, const string& temp, double hden, double metallicity, const Array& radField,
           double ins);

    void perform(CloudyData& data);

private:
    void setup();

    void run();

    void read(CloudyData& data) const;

    // maybe use ions from XRayIonicGasMix and put it here for consistency?
    std::pair<int, int> getElement(string species) const;

    // ================== Input Data ==================
    double _hden;
    double _metallicity;
    Array _radField;
    double _ins;  // mean intensity in W/m2 over the full radiation field range

    // ================== Internal Data ==================
    int _uid;
    string _path;
    const string& _template;
};

#endif
