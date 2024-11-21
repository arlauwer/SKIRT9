#ifndef HEALPIX_HPP
#define HEALPIX_HPP

class HEALPix
{
public:
    static void binHEALPix(double theta, double phi, int& Hi, int& Hj, int Nside)
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
        if (za <= 2. / 3.)
        {
            // equatorial region: all rings have 4*_Nside pixels
            double t1 = Nside * (0.5 + tt);
            double t2 = 0.75 * Nside * z;
            int jp = static_cast<int>(floor(t1 - t2));
            int jm = static_cast<int>(floor(t1 + t2));

            // we first compute the ring index as the offset in the equatorial region (j in [1, 2*_Nside+1])
            Hj = Nside + 1 + jp - jm;
            int kshift = 1 - (Hj & 1);
            int temp = jp + jm + kshift + 1 + 7 * Nside;
            Hi = (Nside > 1) ? ((temp >> 1) & (4 * Nside - 1)) : ((temp >> 1) % (4 * Nside));
            // we now add the correct offset to the ring index to account for the rings in the north polar region
            Hj += Nside - 2;
        }
        else
        {
            // polar region: the number of pixels per ring depends on how far the ring is from the pole
            double tp = tt - static_cast<int>(tt);
            double tmp = (za < 0.99) ? (Nside * sqrt(3. * (1. - za))) : (Nside * sin(theta) / sqrt((1. + za) / 3.));

            int jp = static_cast<int>(tp * tmp);
            int jm = static_cast<int>((1. - tp) * tmp);

            // we first compute the ring index as an offset from the pole (j in [1, _Nside-1])
            Hj = jp + jm + 1;
            Hi = static_cast<int>(tt * Hj);
            // we now convert the ring index to the actual ring index, and distinguish between north and south pole
            if (z < 0)
            {
                Hj = 4 * Nside - Hj - 1;
            }
            else
            {
                --Hj;
            }
        }

        // detect the photon packet in the appropriate pixel of the data cube and at the appropriate distance
        // int l = i + Nx * j;
    }
};

#endif