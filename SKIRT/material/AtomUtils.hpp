/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ATOMUTILS_HPP
#define ATOMUTILS_HPP

#include "Basics.hpp"

/** This static class provides utility functions related to the treatment of atomic and ionic
    species up to atomic number Z=30. */
class AtomUtils final
{
public:
    /** This function returns the atomic number of the specified elements. */
    static int atomToZ(string element);

    /** This function returns the mass of the specified atomic number in SI units. */
    static double mass(int Z);

    /** This function returns the unique index associated with each ion. The index is determined
		from the formula \f$Z(Z+1)/2+N-1\f$. This results in the ions {H+1, H+0, He+2, He+1, ...}
		having indices {0, 1, 2, 3, ...}. */
    static int ionIndex(int Z, int N);

    /** This function returns the full name of an ion given its atomic number and number of
		electrons. e.g. ionName(26, 26) returns "Fe+0". */
    static string ionName(int Z, int N);

    /** This function returns a pair of the atomic number and number of electrons (Z,N) of the
        specified ion string. The ion string can have the following formats, as listed in the table
        below. Any other format will throw an error. <TABLE>
        <TR><TD><B>Format</B></TD><TD><B>Returns</B></TD></TR>
        <TR><TD><TT>'Z'</TT></TD><TD>\f$(Z,N)=(Z,Z)\f$</TD></TR>
        <TR><TD><TT>'Z+'</TT></TD><TD>\f$(Z,N)=(Z,Z-1)\f$</TD></TR>
        <TR><TD><TT>'Z+I'</TT></TD><TD>\f$(Z,N)=(Z,Z-I)\f$</TD></TR> </TABLE> The number of
        electrons, N, can range from 0 to Z, any other value will result in an error. */
    static std::pair<int, int> parseIon(string ion);
};

#endif
