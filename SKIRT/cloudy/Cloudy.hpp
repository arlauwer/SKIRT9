#ifndef CLOUDY_HPP
#define CLOUDY_HPP

#include "Array.hpp"
#include "Constants.hpp"
#include "NR.hpp"
#include <array>
#include <iostream>

namespace cloudy
{
    constexpr int numBins = 24;
    constexpr std::array<double, numBins + 1> edges{
        1.00000000e+04, 6.81292069e+03, 4.64158883e+03, 3.16227766e+03, 2.15443469e+03, 1.46779927e+03, 1.00000000e+03,
        6.81292069e+02, 4.64158883e+02, 3.16227766e+02, 2.15443469e+02, 1.46779927e+02, 1.00000000e+02, 6.81292069e+01,
        4.64158883e+01, 3.16227766e+01, 2.15443469e+01, 1.46779927e+01, 1.00000000e+01, 6.81292069e+00, 4.64158883e+00,
        3.16227766e+00, 2.15443469e+00, 1.46779927e+00, 1.00000000e+00};
    constexpr int numIons = 465;
    constexpr double minRad = 1e-10;
    constexpr double stc = 1e3;                // W/m2 -> erg/s/cm2
    constexpr double Rtm = Constants::iRyd();  // Ryd -> m

    constexpr bool inRange(double Ryd)
    {
        return Ryd >= edges[numBins - 1] && Ryd <= edges[0];
    }
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
