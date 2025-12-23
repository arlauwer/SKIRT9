/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ATOMS_HPP
#define ATOMS_HPP

#include "Constants.hpp"
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
    /** This function returns the atomic number of the specified elements. */
    static inline short atomToZ(string element) { return atomMap.at(element); }

    /** This function returns the atomic number of the specified element. */
    static inline double umass(int Z) { return masses[Z - 1]; }

    static inline double mass(short Z) { return masses[Z - 1] * Constants::amu(); }

    static inline std::pair<short, short> parseIon(string ion)
    {
        // read ions
        ion = StringUtils::squeeze(ion);
        // split ion string into element symbol and ionization number
        std::regex pattern("([A-Za-z]+)\\+?([0-9]*)");
        std::smatch match;

        if (!std::regex_match(ion, match, pattern)) throw FATALERROR("Invalid ion format: " + ion);

        // get ion parameters
        int Z = atomMap.at(match[1].str());
        int N = Z - std::stoi(match[2].str());
        if (N < 1 || N > Z) throw FATALERROR("Invalid ion format: " + ion);
        // add ion parameters
        return std::make_pair(Z, N);
    }
};

#endif
