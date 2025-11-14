#include "Cloudy.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <unordered_map>

////////////////////////////////////////////////////////////////////

namespace
{
    const std::unordered_map<string, int> speciesMap = {
        {"H", 1},   {"He", 2},  {"Li", 3},  {"Be", 4},  {"B", 5},   {"C", 6},   {"N", 7},  {"O", 8},
        {"F", 9},   {"Ne", 10}, {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16},
        {"Cl", 17}, {"Ar", 18}, {"K", 19},  {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24},
        {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30}};
}

////////////////////////////////////////////////////////////////////

CloudyData::CloudyData() : temperature(0.), abundances(0.), opacities(0.), emissivities(0.) {}

////////////////////////////////////////////////////////////////////

Cloudy::Cloudy(int uid, string runsPath, const string& temp, double hden, double metallicity, const Array& radField)
    : _uid(uid), _path(StringUtils::joinPaths(runsPath, StringUtils::toString(uid))), _template(temp), _hden(hden),
      _metallicity(metallicity), _radField(radField)
{}

////////////////////////////////////////////////////////////////////

void Cloudy::perform()
{
    setup();
    run();
    read();
}

////////////////////////////////////////////////////////////////////

void Cloudy::setup()
{
    System::makeDir(_path);  // make uid path

    // probably not needed here
    if (_radField.size() != cloudy::numBins) throw FATALERROR("Cloudy::setup: wrong number of radfield values");

    // sim.in
    string temp = _template;
    temp = StringUtils::replace(temp, "{hden}", StringUtils::toString(_hden));
    temp = StringUtils::replace(temp, "{metallicity}", StringUtils::toString(_metallicity));
    temp = StringUtils::replace(temp, "{ins}", StringUtils::toString(1e-1));
    std::ofstream in = System::ofstream(StringUtils::joinPaths(_path, "sim.in"));
    in.write(temp.data(), temp.size());
    in.close();

    // sed.in
    std::ofstream sed = System::ofstream(StringUtils::joinPaths(_path, "sed.in"));
    // Jnu (~W/m2/Hz) & _radField is (4pi W/m2/m)
    for (int i = 0; i < cloudy::numBins; i++)
    {
        double ryd1 = cloudy::edges[i] * 1.00001;
        double ryd2 = cloudy::edges[i + 1] * 0.9999;
        double lam1 = 9.112662439164599e-08 / ryd1;  // hc / Q / lam / 13.6
        double lam2 = 9.112662439164599e-08 / ryd2;
        double rad = max(_radField[i], 1e-10); // to avoid 0
        double Jnu1 = rad / (4. * M_PI * Constants::c()) * lam1 * lam1;  // same radBin, different lambda
        double Jnu2 = rad / (4. * M_PI * Constants::c()) * lam2 * lam2;  // same radBin, different lambda
        sed << ryd1 << "\t" << Jnu1 << std::endl;                                 // left
        sed << ryd2 << "\t" << Jnu2 << std::endl;                                 // right
    }

    sed.close();
}

////////////////////////////////////////////////////////////////////

void Cloudy::run()
{
    string cmd = "cd " + _path + " && cloudy < sim.in > sim.out";
    System::execute(cmd);
}

////////////////////////////////////////////////////////////////////

void Cloudy::read()
{
    string header;
    string line;

    // temperature
    std::ifstream ovr = System::ifstream(StringUtils::joinPaths(_path, "sim.ovr"));
    getline(ovr, header);  // skip header
    getline(ovr, line);
    ovr.close();
    _data.temperature = StringUtils::toDouble(StringUtils::split(line, "\t")[1]);

    //abundances
    std::ifstream species = System::ifstream(StringUtils::joinPaths(_path, "sim.species"));
    getline(species, header);
    getline(species, line);
    species.close();

    auto speciesHeader = StringUtils::split(header, "\t");
    auto speciesData = StringUtils::split(line, "\t");
    _data.abundances.resize(cloudy::numIons);
    for (int i = 1; i < speciesData.size(); i++)  // skip first column
    {
        auto speciesName = speciesHeader[i];
        int Z, N;
        std::tie(Z, N) = getElement(speciesName);
        if (N > 0)
        {
            int j = Z * (Z - 1) / 2 + (N - 1);
            _data.abundances[j] = StringUtils::toDouble(speciesData[i]);
        }
    }

    // opacities
    std::ifstream opac = System::ifstream(StringUtils::joinPaths(_path, "sim.opac"));
    getline(opac, header);  // skip header
    std::vector<double> opacities;
    for (int i=0; getline(opac, line); i++)
    {
        if (i <  5884 || i > 9024) continue;

        auto cols = StringUtils::split(line, "\t");
        if (cols.size() < 2) continue;
        opacities.push_back(StringUtils::toDouble(cols[2]));
    }
    opac.close();
    _data.opacities = Array(opacities.data(), opacities.size());

    // emissivities
    std::ifstream emis = System::ifstream(StringUtils::joinPaths(_path, "sim.con"));
    getline(emis, header);  // skip header
    std::vector<double> emissivities;
    for (int i=0; getline(emis, line); i++)
    {
        if (i <  5884 || i > 9024) continue;

        auto cols = StringUtils::split(line, "\t");
        if (cols.size() < 3) continue;
        // 4pi erg/s/cm2
        emissivities.push_back(StringUtils::toDouble(cols[3]));
    }
    emis.close();
    _data.emissivities = Array(emissivities.data(), emissivities.size());
}

////////////////////////////////////////////////////////////////////

std::pair<int, int> Cloudy::getElement(string species)
{
    auto parts = StringUtils::split(species, "+");

    if (speciesMap.find(parts[0]) == speciesMap.end()) return {0, 0};

    int Z = speciesMap.at(parts[0]);
    int N = Z;
    if (parts.size() == 2)
    {
        if (parts[1].size() == 0)
            N = Z - 1;
        else
            N = Z - StringUtils::toInt(parts[1]);
    }
    return {Z, N};
}

////////////////////////////////////////////////////////////////////
