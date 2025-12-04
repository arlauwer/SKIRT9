#include "Cloudy.hpp"
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

CloudyData::CloudyData() : temperature(0.), abundances(0., 465), opacities(0., 3141), emissivities(0., 3141) {}

////////////////////////////////////////////////////////////////////

Cloudy::Cloudy(int uid, string runsPath, const string& temp, double hden, double metallicity, const Array& radField,
               double ins)
    : _hden(hden * 1e-6), _metallicity(metallicity), _radField(radField), _ins(ins * cloudy::stc), _uid(uid),
      _path(StringUtils::joinPaths(runsPath, StringUtils::toString(uid))), _template(temp)
{
    // hden (1/m3) -> _hden (1/cm3)
    // metallicity (1) -> _metallicity (1)
    // rad (W/4pi/m2/m) -> _rad (~W/m2/m)
    // ins (W/m2) -> _ins (erg/s/cm2)
}

////////////////////////////////////////////////////////////////////

void Cloudy::perform(CloudyData& data)
{
    setup();
    run();
    read(data);
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
    temp = StringUtils::replace(temp, "{ins}", StringUtils::toString(_ins));
    std::ofstream in = System::ofstream(StringUtils::joinPaths(_path, "sim.in"));
    in.write(temp.data(), temp.size());
    in.close();

    // sed.in
    std::ofstream sed = System::ofstream(StringUtils::joinPaths(_path, "sed.in"));
    // Jlambda (~W/m2/m) & _radField is (4pi W/m2/m)
    for (int i = 0; i < cloudy::numBins; i++)
    {
        double left = cloudy::edges[i] * 0.9999;
        double right = cloudy::edges[i + 1] * 1.00001;
        double rad = max(_radField[i], cloudy::minRad);  // to avoid 0

        if (i == 0)
        {
            sed << left << "\t" << rad << " Flambda\n";
            sed << right << "\t" << rad << "\n";
        }
        else
        {
            sed << left << "\t" << rad << "\n";
            sed << right << "\t" << rad << "\n";
        }
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

void Cloudy::read(CloudyData& data) const
{
    string header;
    string line;

    // temperature
    std::ifstream ovr = System::ifstream(StringUtils::joinPaths(_path, "sim.ovr"));
    getline(ovr, header);  // skip header
    getline(ovr, line);
    ovr.close();
    data.temperature = StringUtils::toDouble(StringUtils::split(line, "\t")[1]);

    //abundances
    std::ifstream species = System::ifstream(StringUtils::joinPaths(_path, "sim.species"));
    getline(species, header);
    getline(species, line);
    species.close();

    auto speciesHeader = StringUtils::split(header, "\t");
    auto speciesData = StringUtils::split(line, "\t");
    data.abundances.resize(cloudy::numIons);
    for (size_t i = 1; i < speciesData.size(); i++)  // skip first column
    {
        auto speciesName = speciesHeader[i];
        int Z, N;
        std::tie(Z, N) = getElement(speciesName);
        if (N > 0)
        {
            int j = Z * (Z - 1) / 2 + (N - 1);
            data.abundances[j] = StringUtils::toDouble(speciesData[i]) * 1e6;  // 1/cm3 -> 1/m3
        }
    }

    // opacities
    std::ifstream opac = System::ifstream(StringUtils::joinPaths(_path, "sim.opac"));
    getline(opac, header);  // skip header
    std::vector<double> opacities;
    for (int i = 0; getline(opac, line); i++)
    {
        auto cols = StringUtils::split(line, "\t");
        if (cols.size() < 2) continue;

        if (i < 5884 || i > 9024) continue;

        opacities.push_back(StringUtils::toDouble(cols[2]));
    }
    opac.close();
    std::reverse(opacities.begin(), opacities.end());
    data.opacities = Array(opacities.data(), opacities.size());
    if (data.opacities.size() != 3141) throw FATALERROR("Cloudy::read: wrong number of opacities");

    // emissivities
    std::ifstream emis = System::ifstream(StringUtils::joinPaths(_path, "sim.con"));
    getline(emis, header);  // skip header
    std::vector<double> emissivities;
    for (int i = 0; getline(emis, line); i++)
    {
        auto cols = StringUtils::split(line, "\t");
        if (cols.size() < 3) continue;

        if (i < 5884 || i > 9024) continue;

        // 4pi erg/s/cm2
        // should be W/m2/m
        emissivities.push_back(StringUtils::toDouble(cols[3]));
    }
    emis.close();
    std::reverse(emissivities.begin(), emissivities.end());
    data.emissivities = Array(emissivities.data(), emissivities.size());
    if (data.emissivities.size() != 3141) throw FATALERROR("Cloudy::read: wrong number of emissivities");
}

////////////////////////////////////////////////////////////////////

std::pair<int, int> Cloudy::getElement(string species) const
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
