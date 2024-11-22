#include "HEALPix.hpp"
#include "math.h"

template<typename I> inline I imodulo(I v1, I v2)
{
    return (v1 >= 0) ? ((v1 < v2) ? v1 : (v1 % v2)) : ((v1 % v2) + v2);
}

HEALPix::HEALPix() {}

HEALPix::HEALPix(int order)
    : _order(order), _Nside(1 << _order), _Nbins(12 * _Nside << _order), _Ncap((_Nside << _order - _Nside) << 1)
{}

int HEALPix::binHEALPix(double theta, double phi) const
{
    // p.spherical() returns theta in [0, pi] and phi in [-pi, pi]
    // the HEALPix mapping algorithm expects theta in [0, pi] and phi in [0, 2*pi]
    // we center the output image on the crosshair direction by offsetting phi
    phi += M_PI;

    // the code below was mostly copied from healpix_base.cc in the official HEALPix repository
    // it corresponds to the loc2pix function using the RING scheme
    // unlike the loc2pix function, we do not compute a full RING pixel index, but restrict
    // ourselves to a computation of the ring index and pixel-in-ring index, which we respectively
    // use as vertical and horizontal pixel index
    double z = cos(theta);
    double za = abs(z);
    double tt = fmod(2. * phi / M_PI, 4.);

    // i: pixel-in-ring index (horizontal pixel index)
    // j: ring index (vertical pixel index)
    // first figure out if we are in the equatorial or polar region
    if (za <= 2. / 3.)  // Equatorial region
    {
        double temp1 = _Nside * (0.5 + tt);
        double temp2 = _Nside * z * 0.75;
        int jp = int(temp1 - temp2);  // index of  ascending edge line
        int jm = int(temp1 + temp2);  // index of descending edge line

        // ring number counted from z=2/3
        int ir = _Nside + 1 + jp - jm;  // in {1,2n+1}
        int kshift = 1 - (ir & 1);      // kshift=1 if ir even, 0 otherwise

        int ip = (jp + jm - _Nside + kshift + 1) / 2;  // in {0,4n-1}
        ip = imodulo(ip, 4 * _Nside);

        return _Ncap + (ir - 1) * 4 * _Nside + ip;
    }
    else  // North & South polar caps
    {
        double tp = tt - int(tt);
        double tmp = _Nside * sqrt(3 * (1 - za));

        int jp = int(tp * tmp);          // increasing edge line index
        int jm = int((1.0 - tp) * tmp);  // decreasing edge line index

        int ir = jp + jm + 1;   // ring number counted from the closest pole
        int ip = int(tt * ir);  // in {0,4*ir-1}
        ip = imodulo(ip, 4 * ir);

        if (z > 0)
            return 2 * ir * (ir - 1) + ip;
        else
            return _Nbins - 2 * ir * (ir + 1) + ip;
    }
}
