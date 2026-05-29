/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ATOMS_HPP
#define ATOMS_HPP

#include "FatalError.hpp"
#include "StringUtils.hpp"
#include <map>
#include <regex>

////////////////////////////////////////////////////////////////////

namespace
{
    static const std::map<string, short> atomMap = {
        {"H", 1},   {"He", 2},  {"Li", 3},  {"Be", 4},  {"B", 5},   {"C", 6},   {"N", 7},  {"O", 8},
        {"F", 9},   {"Ne", 10}, {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16},
        {"Cl", 17}, {"Ar", 18}, {"K", 19},  {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24},
        {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30}};

    static const vector<double> masses = {1.0079,  4.0026, 6.941,   9.01218, 10.81,   12.011,  14.0067, 15.9994,
                                          18.9984, 20.179, 22.9898, 24.305,  26.9815, 28.0855, 30.9738, 32.06,
                                          35.453,  39.948, 39.0983, 40.08,   44.9559, 47.9,    50.9415, 51.996,
                                          54.938,  55.847, 58.9332, 58.7,    63.546,  65.38};
}

/** This class provides utility functions related to the treatment of atomic and ionic species. */
class Atoms
{
public:
    struct Ion
    {
        int Z;  // atomic number
        int N;  // number of electrons
    };

    // Ions are ordered as such
    // i    =   0,    1,    2,    3,    4,    5, ...
    // name = H+1,  H+0, He+2, He+1, He+0, Li+3, ...
    // Z-N  = 1-0, 1-1,   2-0,  2-1,  2-2,  3-0, ...
    static inline int ionIndex(int Z, int N) { return Z * (Z + 1) / 2 + N - 1; }
    static inline int ionIndex(Ion ion) { return ionIndex(ion.Z, ion.N); }

    template<int numIons, int numAtoms> constexpr static std::array<Ion, numIons> initIons()
    {
        std::array<Ion, numIons> ions{};
        for (int Z = 1; Z <= numAtoms; Z++)
        {
            for (int N = 0; N <= Z; N++)
            {
                int i = ionIndex(Z, N);
                ions[i].Z = Z;
                ions[i].N = N;
            }
        }
        return ions;
    }

    static string ionName(int Z, int N)
    {
        auto it = std::find_if(atomMap.begin(), atomMap.end(), [Z](const auto& pair) { return pair.second == Z; });
        if (it == atomMap.end()) throw FATALERROR("Could not find element for Z = " + std::to_string(Z));

        if (N < 0 || N > Z)
            throw FATALERROR("Invalid ionization number for (Z,N) = (" + std::to_string(Z) + "," + std::to_string(N)
                             + ")");

        return it->first + "+" + std::to_string(Z - N);
    }

    /** This function returns the atomic number of the specified elements. */
    static inline short atomToZ(string element) { return atomMap.at(element); }

    /** This function returns the atomic number of the specified element. */
    static inline double mass(int Z) { return masses[Z - 1]; }

    static inline Ion parseIon(string ion)
    {
        // read ions
        ion = StringUtils::squeeze(ion);
        // split ion string into element symbol and ionization number
        std::regex pattern("([A-Za-z]+)(\\+?)([0-9]*)");
        std::smatch match;

        if (!std::regex_match(ion, match, pattern)) throw FATALERROR("Could not parse ion format: " + ion);

        int Z, N;
        string element = match[1].str();
        string plus = match[2].str();
        string number = match[3].str();

        Z = atomMap.at(element);

        if (number.empty() && plus.empty())
            N = Z;
        else if (number.empty())
            N = Z - 1;
        else
            N = Z - StringUtils::toInt(number);

        return Ion{Z, N};
    }
};

#endif
