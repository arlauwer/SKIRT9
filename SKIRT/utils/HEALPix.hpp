#ifndef HEALPIX_HPP
#define HEALPIX_HPP

/**
 * HEALPix binning class that also keeps information about how many bins. If order is set to -1 there is just 1 bin which is the entire sphere.
 */
class HEALPix
{
public:
    HEALPix();

    HEALPix(int order);

    int binHEALPix(double theta, double phi) const;

    int Nbins() const { return _Nbins; }

private:
    int _order;
    int _Nside;
    int _Nbins{1};  // default value = 1 bin for all sky
    int _Ncap;
};

#endif