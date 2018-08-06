/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SCATTERINGMIXINTERFACE_HPP
#define SCATTERINGMIXINTERFACE_HPP

#include "Direction.hpp"
class PhotonPacket;
class StokesVector;

////////////////////////////////////////////////////////////////////

/** The ScatteringMixInterface interface offers functions to handle scattering events involving an
    incoming photon packet and its outgoing direction. When polarization is not supported (as
    indicated by the hasPolarization() function), these scattering functions may internally employ
    approximations such as using the Henyey-Greenstein function parameterized by an asymmetry
    parameter \f$g_\lambda\f$. When polarization is supported, these functions should obtain the
    phase function from the Mueller matrix in combination with calculating the Stokes vector
    updates. */
class ScatteringMixInterface
{
    //======== Construction - Destruction =======

protected:
    /** The empty constructor for the interface. */
    ScatteringMixInterface() { }

public:
    /** The empty destructor for the interface. */
    virtual ~ScatteringMixInterface() { }

    //======== Functions for handling scattering events =======

public:
    /** This function returns true if this material mix supports polarization by scattering; false
        otherwise. */
    virtual bool hasScatteringPolarization() const = 0;

    /** This function generates a new direction \f${\bf{k}}_{\text{new}}\f$ in case the specified
        photon packet scatters, and calculates the new polarization state of the scattered photon
        packet. The function passes the new direction to the caller as its return value, and stores
        the new polarization state in the provided Stokes vector. It is permitted for the provided
        Stokes vector to actually reside in the specified photon packet. For a material mix that
        doesn't support polarization, the provided Stokes vector is not modified. */
    virtual Direction scatteringDirectionAndPolarization(StokesVector* out, const PhotonPacket* pp) const = 0;

    /** This function calculates the polarization state appropriate for a peel off photon packet
        generated by a scattering event for the specified photon packet, and stores the result in
        the provided Stokes vector. The resulting Stokes vector must be rotated so that it is
        properly oriented for the reference frame of the instrument to which the peel-off photon
        packet is headed. For a material mix that doesn't support polarization, the provided Stokes
        vector is not modified. */
    virtual void scatteringPeelOffPolarization(StokesVector* out, const PhotonPacket* pp, Direction bfknew,
                                               Direction bfkx, Direction bfky) = 0;

    /** This function returns the value of the scattering phase function in case the specified
        photon packet is scattered to the specified new direction, where the phase function is
        normalized as \f[\int\Phi_\lambda(\Omega)\,\mathrm{d}\Omega=4\pi.\f] */
    virtual double phaseFunctionValue(const PhotonPacket* pp, Direction bfknew) const = 0;
};

////////////////////////////////////////////////////////////////////

#endif