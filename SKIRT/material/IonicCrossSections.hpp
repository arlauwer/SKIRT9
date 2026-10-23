/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef IONICCROSSSECTIONS_HPP
#define IONICCROSSSECTIONS_HPP

#include "AtomUtils.hpp"
#include "Constants.hpp"
#include "Range.hpp"
#include "TextInFile.hpp"
#include <cmath>
#include <unordered_set>

////////////////////////////////////////////////////////////////////

/** The IonicCrossSections namespace offers the photo-absorption and fluorescence data for
    individual ions that is shared by the XRayIonicGasMix and XRayCloudyGasMix classes, and the
    functions to load this data from the resource files.

    The photo-absorption data, loaded from the resource file <tt>Ionic_PA.txt</tt>, holds the
    analytic fits of Verner & Yakovlev 1995 (A&AS, 109, 125) for the partial photo-ionization cross
    section of each subshell of each ion with \f$1\le Z\le 30\f$. The fluorescence data, loaded
    from the resource file <tt>Ionic_FL.txt</tt>, holds for each subshell of each ion the
    fluorescence lines that may follow the creation of a vacancy in that subshell, each with its
    yield, its central energy and, for some lines, the width of its Lorentz profile. The yields
    include the photons emitted during the subsequent cascade, so that the total yield for a
    subshell may exceed unity.

    The cross sections are nonzero only between 4.3 eV and 500 keV (see nonZeroRange). */
namespace IonicCrossSections
{
    /** This function converts a photon energy in eV to a wavelength in m, and vice versa (the
        conversion is the same in both directions). */
    constexpr double wavelengthToFromEnergy(double x)
    {
        constexpr double front = Constants::h() * Constants::c() / Constants::Qelectron();
        return front / x;
    }

    /** The wavelength range over which the cross sections may be nonzero. */
    constexpr Range nonZeroRange(wavelengthToFromEnergy(500e3), wavelengthToFromEnergy(4.3));

    /** This structure holds the photo-absorption fit parameters for a single subshell of a single
        ion, and offers functions to evaluate the corresponding cross section. */
    struct PhotoAbsorbResource
    {
        PhotoAbsorbResource(const Array& a)
            : Z(a[0]), N(a[1]), n(a[2]), l(a[3]), Eth(a[4]), E0(a[5]), sigma0(a[6]), ya(a[7]), P(a[8]), yw(a[9])
        {}

        int ionIndex{-1};                    // index of the ion in the client's ion list
        int Z;                               // atomic number
        int N;                               // number of electrons
        int n;                               // principal quantum number of the shell
        int l;                               // orbital quantum number of the subshell
        double Eth;                          // subshell ionization threshold energy (eV)
        double E0;                           // fit parameter (eV)
        double sigma0;                       // fit parameter (Mb = 10^-22 m^2)
        double ya;                           // fit parameter (1)
        double P;                            // fit parameter (1)
        double yw;                           // fit parameter (1)
        static constexpr double Emax = 5e5;  // maximum energy for validity of the formula (eV)

        double Es{0.};        // thermal energy dispersion at the threshold energy (eV); set by the client
        double sigmamax{0.};  // cross section at the threshold energy plus 2*Es (m2); set by the client

        /** This function returns the photo-absorption cross section in m2 for the given energy in
            eV, without taking into account thermal dispersion. */
        double photoAbsorbSection(double E) const
        {
            if (E < Eth || E >= Emax) return 0.;

            double x = E / E0;
            double y = x;
            double xm1 = x - 1.;
            double Q = 5.5 + l - 0.5 * P;
            double F = (xm1 * xm1 + yw * yw) * std::pow(y, -Q) * std::pow(1. + std::sqrt(y / ya), -P);
            return 1e-22 * sigma0 * F;  // from Mb to m2
        }

        /** This function returns the photo-absorption cross section in m2 for the given energy in
            eV, approximating thermal dispersion by replacing the steep threshold transition by a
            sigmoid error function with the parameters \em Es and \em sigmamax, which must have
            been set by the client. */
        double photoAbsorbThermalSection(double E) const
        {
            if (E <= Eth - 2. * Es) return 0.;
            if (E >= Eth + 2. * Es) return photoAbsorbSection(E);
            return sigmamax * (0.5 + 0.5 * std::erf((E - Eth) / Es));
        }
    };

    /** This structure holds the parameters of a single fluorescence line that may follow the
        creation of a vacancy in a given subshell of a given ion. */
    struct FluorescenceResource
    {
        FluorescenceResource(const Array& a) : Z(a[0]), N(a[1]), n(a[2]), l(a[3]), omega(a[4]), E(a[5]), W(a[6]) {}

        int paIndex{-1};  // index of the photo-absorption transition; set by the client
        int Z;            // atomic number
        int N;            // number of electrons
        int n;            // principal quantum number of the shell with the hole
        int l;            // orbital quantum number of the subshell with the hole
        double omega;     // fluorescence yield (1)
        double E;         // (central) energy of the emitted photon (eV)
        double W;         // FWHM of the Lorentz shape for the emitted photon (eV), or zero
    };

    /** This function loads data from a resource file, only keeping the data for the ions in the
        specified set. The data is loaded into a vector of structures of type S that can be
        constructed from an array with C elements, i.e. C columns. The first two columns of the
        data must always be the ion parameters Z and N. The set holds the indices of the ions to
        be kept, as calculated by the AtomUtils::ionIndex() function. */
    template<class S, int C>
    vector<S> loadPresent(SimulationItem* item, string filename, string description,
                          const std::unordered_set<int>& ionSet)
    {
        vector<S> result;
        TextInFile infile(item, filename, description, true);
        infile.addColumn("Z");
        infile.addColumn("N");
        for (int i = 2; i < C; i++) infile.addColumn(string());

        Array row;
        while (infile.readRow(row))
        {
            int Z = row[0];
            int N = row[1];
            int index = AtomUtils::ionIndex(Z, N);

            // only keep ions that are present
            if (ionSet.count(index)) result.emplace_back(row);
        }
        return result;
    }

    /** This function loads the photo-absorption data for the ions in the specified set from the
        resource file <tt>Ionic_PA.txt</tt>, in the order of that file. */
    vector<PhotoAbsorbResource> loadPhotoAbsorption(SimulationItem* item, const std::unordered_set<int>& ionSet);

    /** This function loads the fluorescence data for the ions in the specified set from the
        resource file <tt>Ionic_FL.txt</tt>, in the order of that file. */
    vector<FluorescenceResource> loadFluorescence(SimulationItem* item, const std::unordered_set<int>& ionSet);
}

////////////////////////////////////////////////////////////////////

#endif
