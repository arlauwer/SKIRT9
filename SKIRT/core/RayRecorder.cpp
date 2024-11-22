/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "RayRecorder.hpp"
#include "FITSInOut.hpp"
#include "Indices.hpp"
#include "LockFree.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PathSegmentGenerator.hpp"
#include "PhotonPacket.hpp"
#include "ProcessManager.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // indices for detector arrays that need calibration
    enum {
        Total = 0,
        Transparent,
        PrimaryDirect,
        PrimaryScattered,
        SecondaryTransparent,
        SecondaryDirect,
        SecondaryScattered,
        TotalQ,
        TotalU,
        TotalV,
        PrimaryScatteredLevel
    };

    // the highest contribution power to track, i.e. the largest k in sum(w_i**k)
    //  - include k=0, which tracks the number of detections
    //  - thus, the number of detector arrays for statistics is this number plus one
    //  - these detector arrays do not need calibration!
    const int maxContributionPower = 4;
}

////////////////////////////////////////////////////////////////////

RayRecorder::RayRecorder(const SimulationItem* parentItem) : _parentItem{parentItem} {}

////////////////////////////////////////////////////////////////////

void RayRecorder::setSimulationInfo(string instrumentName, const WavelengthGrid* lambdagrid, bool hasMedium,
                                    bool hasMediumEmission)
{
    _instrumentName = instrumentName;
    _lambdagrid = lambdagrid;
    _hasMedium = hasMedium;
    _hasMediumEmission = hasMediumEmission;
}

////////////////////////////////////////////////////////////////////

void RayRecorder::setUserFlags(bool recordComponents, int numScatteringLevels, bool recordPolarization,
                               bool recordStatistics)
{
    _recordComponents = recordComponents;
    _numScatteringLevels = numScatteringLevels;
    _recordPolarization = recordPolarization;
    _recordStatistics = recordStatistics;
}

////////////////////////////////////////////////////////////////////

void RayRecorder::setObserverAngles(double inclination, double azimuth, double roll)
{
    _inclination = inclination;
    _azimuth = azimuth;
    _roll = roll;
}

////////////////////////////////////////////////////////////////////

void RayRecorder::setRestFrameDistance(double distance)
{
    _redshift = 0.;
    _angularDiameterDistance = distance;
    _luminosityDistance = distance;
}

////////////////////////////////////////////////////////////////////

void RayRecorder::setObserverFrameRedshift(double redshift, double angularDiameterDistance, double luminosityDistance)
{
    _redshift = redshift;
    _angularDiameterDistance = angularDiameterDistance;
    _luminosityDistance = luminosityDistance;
}

////////////////////////////////////////////////////////////////////

void RayRecorder::includeFluxDensityForDistant()
{
    _includeFluxDensity = true;
}

////////////////////////////////////////////////////////////////////

void RayRecorder::includeSurfaceBrightnessForDistant(int numPixelsX, int numPixelsY, double pixelSizeX,
                                                     double pixelSizeY, double centerX, double centerY)
{
    _includeSurfaceBrightness = true;
    _numPixelsX = numPixelsX;
    _numPixelsY = numPixelsY;
    _pixelSizeX = pixelSizeX;
    _pixelSizeY = pixelSizeY;
    _centerX = centerX;
    _centerY = centerY;
}

////////////////////////////////////////////////////////////////////

void RayRecorder::includeSurfaceBrightnessForLocal(int numPixelsX, int numPixelsY, double solidAnglePerPixel,
                                                   double incrementX, double incrementY, double centerX, double centerY,
                                                   string quantityXY)
{
    _includeSurfaceBrightness = true;
    _local = true;
    _numPixelsX = numPixelsX;
    _numPixelsY = numPixelsY;
    _solidAnglePerPixel = solidAnglePerPixel;
    _incrementX = incrementX;
    _incrementY = incrementY;
    _centerX = centerX;
    _centerY = centerY;
    _quantityXY = quantityXY;
}

////////////////////////////////////////////////////////////////////

void RayRecorder::finalizeConfiguration()
{
    // get a pointer to the medium system, if present
    _ms = _parentItem->find<MediumSystem>(false);

    // get array lengths
    _numPixelsInFrame = _numPixelsX * _numPixelsY;  // convert to size_t before calculating lenIFU
    size_t lenSED = _includeFluxDensity ? _lambdagrid->numBins() : 0;
    size_t lenIFU = _includeSurfaceBrightness ? _numPixelsInFrame * _lambdagrid->numBins() : 0;

    // do not try to record components if there is no medium
    _recordTotalOnly = !_recordComponents || !_hasMedium;

    // allocate the appropriate number of flux detector arrays
    _sed.resize(PrimaryScatteredLevel + _numScatteringLevels);
    _ifu.resize(PrimaryScatteredLevel + _numScatteringLevels);

    // resize the flux detector arrays according to the configuration
    if (_recordTotalOnly)
    {
        _sed[Total].resize(lenSED);
        _ifu[Total].resize(lenIFU);
    }
    else
    {
        _sed[Transparent].resize(lenSED);
        _ifu[Transparent].resize(lenIFU);
        _sed[PrimaryDirect].resize(lenSED);
        _ifu[PrimaryDirect].resize(lenIFU);
        _sed[PrimaryScattered].resize(lenSED);
        _ifu[PrimaryScattered].resize(lenIFU);

        for (int i = 0; i != _numScatteringLevels; ++i)
        {
            _sed[PrimaryScatteredLevel + i].resize(lenSED);
            _ifu[PrimaryScatteredLevel + i].resize(lenIFU);
        }
        if (_hasMediumEmission)
        {
            _sed[SecondaryTransparent].resize(lenSED);
            _ifu[SecondaryTransparent].resize(lenIFU);
            _sed[SecondaryDirect].resize(lenSED);
            _ifu[SecondaryDirect].resize(lenIFU);
            _sed[SecondaryScattered].resize(lenSED);
            _ifu[SecondaryScattered].resize(lenIFU);
        }
    }
    if (_recordPolarization)
    {
        _sed[TotalQ].resize(lenSED);
        _ifu[TotalQ].resize(lenIFU);
        _sed[TotalU].resize(lenSED);
        _ifu[TotalU].resize(lenIFU);
        _sed[TotalV].resize(lenSED);
        _ifu[TotalV].resize(lenIFU);
    }

    // allocate and resize the statistics detector arrays
    if (_recordStatistics)
    {
        _wsed.resize(maxContributionPower + 1);
        _wifu.resize(maxContributionPower + 1);
        for (auto& array : _wsed) array.resize(lenSED);
        for (auto& array : _wifu) array.resize(lenIFU);
    }

    // calculate and log allocated memory size
    size_t allocatedSize = 0;
    for (const auto& array : _sed) allocatedSize += array.size();
    for (const auto& array : _ifu) allocatedSize += array.size();
    for (const auto& array : _wsed) allocatedSize += array.size();
    for (const auto& array : _wifu) allocatedSize += array.size();
    _parentItem->find<Log>()->info(_parentItem->typeAndName() + " allocated "
                                   + StringUtils::toMemSizeString(allocatedSize * sizeof(double)) + " of memory");
}

void RayRecorder::rayTrace()
{
    // get frame configuration
    int Nxp = _numPixelsX;
    int Nyp = _numPixelsY;
    double xcent = _centerX;
    double ycent = _centerY;
    double xfov = _numPixelsX * _pixelSizeX;
    double yfov = _numPixelsY * _pixelSizeY;
    double xpmin = _centerX - 0.5 * xfov;
    double xpsiz = _pixelSizeX;
    double ypmin = _centerY - 0.5 * yfov;
    double ypsiz = _pixelSizeY;

    // get the oversampling rate
    int Nsampling = 1;

    // calculate sines and cosines for the transformation from observer to model coordinates
    double costheta = cos(_inclination);
    double sintheta = sin(_inclination);
    double cosphi = cos(_azimuth);
    double sinphi = sin(_azimuth);
    double cosomega = cos(_roll);
    double sinomega = sin(_roll);

    // determine the observer axes as directions in model coordinates (inverse transform of frame instrument)
    // k_z is the direction from observer to model (i.e. left-handed coordinate frame)
    Direction kz(-cosphi * sintheta, -sinphi * sintheta, -costheta, false);
    Direction ky(-cosphi * costheta * cosomega - sinphi * sinomega, -sinphi * costheta * cosomega + cosphi * sinomega,
                 sintheta * cosomega, false);
    Direction kx(cosphi * costheta * sinomega - sinphi * cosomega, sinphi * costheta * sinomega + cosphi * cosomega,
                 -sintheta * sinomega, false);

    // get a distance that should be well outside of the model (but with a similar order of magnitude)
    double zp = -10. * (xfov + yfov); // negative sign

    // allocate result array with the appropriate size and initialize contents to zero
    int numWavelengths = _lambdagrid->numBins();
    Table<3> vvv(Nyp, Nxp, numWavelengths);

    // calculate the results in parallel
    auto log = _parentItem->find<Log>();
    auto ms = _parentItem->find<MediumSystem>();
    auto grid = ms->grid();
    log->infoSetElapsed(Nyp);
    auto parallel = _parentItem->find<ParallelFactory>()->parallelDistributed();

    // find HEALPix bin
    kz = -kz;
    double theta, phi;
    kz.spherical(theta, phi);
    // theta += M_PI;  // done by kz = -kz
    int Hi = ms->radiationFieldDirectionBin().binHEALPix(theta, phi);
    // ms->test();

    parallel->call(Nyp, [&vvv, &ms, &grid, xpmin, xpsiz, ypmin, ypsiz, kx, ky, kz, zp, costheta, sintheta, cosphi,
                         sinphi, cosomega, sinomega, Nxp, Hi, log, Nsampling,
                         numWavelengths](size_t firstIndex, size_t numIndices) {
        int Nsampling2 = Nsampling * Nsampling;

        // loop over pixels
        for (size_t j = firstIndex; j != firstIndex + numIndices; ++j)
        {
            for (int i = 0; i != Nxp; ++i)
            {
                for (int is = 0; is < Nsampling; ++is)
                {
                    for (int js = 0; js < Nsampling; ++js)
                    {
                        // transform pixel indices to observer coordinates
                        double xp = xpmin + (i + (is + 1.0) / (Nsampling + 1.0)) * xpsiz;
                        double yp = ypmin + (j + (js + 1.0) / (Nsampling + 1.0)) * ypsiz;

                        // transform observer coordinates to model coordinates (inverse transform of frame instrument)
                        double xpp = sinomega * xp - cosomega * yp;
                        double ypp = cosomega * xp + sinomega * yp;
                        double zpp = zp;
                        double x = cosphi * costheta * xpp - sinphi * ypp + cosphi * sintheta * zpp;
                        double y = sinphi * costheta * xpp + cosphi * ypp + sinphi * sintheta * zpp;
                        double z = -sintheta * xpp + costheta * zpp;

                        // loop over wavelengths
                        for (int ell = 0; ell < numWavelengths; ++ell)
                        {
                            // get a segment generator and initialize the path
                            auto generator = grid->createPathSegmentGenerator();
                            SpatialGridPath path(Position(x, y, z), kz);
                            generator->start(&path);

                            double I = 0.;  // accumulated intensity

                            // accumulate values along the path
                            while (generator->next())
                            {
                                int m = generator->m();
                                if (m < 0) continue;

                                Position pppse = grid->centralPositionInCell(m);

                                double j =
                                    ms->specificIntensity(m, Hi)[ell] / generator->ds();  // volume emissivity
                                // solve RT equation for constant emissivity and opacity
                                double k = ms->opacityExt(ell, m);
                                double ext = exp(-k * generator->ds());

                                if (k)
                                {
                                    I = I * ext + j / k * (1. - ext);
                                }
                                else
                                {
                                    I += j;
                                }
                            }
                            vvv(j, i, ell) += I / Nsampling2;
                        }
                    }
                }
            }
        }
    });
    ProcessManager::sumToRoot(vvv.data(), true);

    // Loop over pixels and wavelengths
    for (size_t j = 0; j < Nyp; ++j)
    {
        for (int i = 0; i < Nxp; ++i)
        {
            for (int ell = 0; ell < numWavelengths; ++ell)
            {
                double I = vvv(j, i, ell);

                // Record in IFU arrays
                if (_includeSurfaceBrightness)
                {
                    size_t lell = j * Nxp + i + ell * _numPixelsInFrame;

                    if (_recordTotalOnly)
                    {
                        LockFree::add(_ifu[Total][lell], I);
                    }
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

void RayRecorder::calibrateAndWrite()
{
    // collect recorded data from all processes
    for (auto& array : _sed) ProcessManager::sumToRoot(array);
    for (auto& array : _ifu) ProcessManager::sumToRoot(array);
    for (auto& array : _wsed) ProcessManager::sumToRoot(array);
    for (auto& array : _wifu) ProcessManager::sumToRoot(array);

    // calibrate and write only in the root process
    if (!ProcessManager::isRoot()) return;

    // calculate front factors for converting from recorded quantities to output quantities
    // (for local instruments, the distance correction already happened)
    double fourpid2 = 4. * M_PI * (_local ? 1. : _luminosityDistance * _luminosityDistance);
    double omega = _local ? _solidAnglePerPixel
                          : 4. * atan(0.5 * _pixelSizeX / _angularDiameterDistance)
                                * atan(0.5 * _pixelSizeY / _angularDiameterDistance);

    // convert from recorded quantities to output quantities and from internal units to user-selected output units
    // (for performance reasons, determine the units scaling factor only once for each wavelength)
    Units* units = _parentItem->find<Units>();
    int numWavelengths = _lambdagrid->numBins();
    for (int ell = 0; ell != numWavelengths; ++ell)
    {
        // SEDs
        if (_includeFluxDensity)
        {
            double factor = 1. / fourpid2 / _lambdagrid->effectiveWidth(ell)
                            * units->ofluxdensity(_lambdagrid->wavelength(ell), 1.);
            for (auto& array : _sed)
                if (array.size()) array[ell] *= factor;
        }
        // IFUs
        if (_includeSurfaceBrightness)
        {
            double factor = 1. / fourpid2 / omega / _lambdagrid->effectiveWidth(ell)
                            * units->osurfacebrightness(_lambdagrid->wavelength(ell), 1.);
            size_t begin = ell * _numPixelsInFrame;
            size_t end = begin + _numPixelsInFrame;
            for (auto& array : _ifu)
                if (array.size())
                    for (size_t lell = begin; lell != end; ++lell) array[lell] *= factor;
        }
    }

    // write SEDs to a single text file (with multiple columns)
    if (_includeFluxDensity)
    {
        // Build a list of column names and corresponding pointers to sed arrays (which may be empty)
        vector<string> sedNames;
        vector<Array*> sedArrays;

        // add the total flux; if we didn't record it directly, calculate it now
        sedNames.push_back("total flux");
        Array sedTotal;
        if (_recordTotalOnly)
            sedArrays.push_back(&_sed[Total]);
        else
        {
            sedTotal = _sed[PrimaryDirect] + _sed[PrimaryScattered];
            if (_hasMediumEmission) sedTotal += _sed[SecondaryDirect] + _sed[SecondaryScattered];
            sedArrays.push_back(&sedTotal);
        }

        // add the flux components, if requested
        // we always add all of them, even if some of them are zero
        if (_recordComponents)
        {
            // add transparent flux
            // if we did not actually record components (because there are no media), use the total flux instead
            sedNames.push_back("transparent flux");
            sedArrays.push_back(_recordTotalOnly ? &_sed[Total] : &_sed[Transparent]);

            // add the actual components of the total flux
            sedNames.insert(sedNames.end(), {"direct primary flux", "scattered primary flux", "direct secondary flux",
                                             "scattered secondary flux", "transparent secondary flux"});
            sedArrays.insert(sedArrays.end(), {&_sed[PrimaryDirect], &_sed[PrimaryScattered], &_sed[SecondaryDirect],
                                               &_sed[SecondaryScattered], &_sed[SecondaryTransparent]});
        }

        // add the polarization components, if requested
        if (_recordPolarization)
        {
            sedNames.insert(sedNames.end(), {"total Stokes Q", "total Stokes U", "total Stokes V"});
            sedArrays.insert(sedArrays.end(), {&_sed[TotalQ], &_sed[TotalU], &_sed[TotalV]});
        }

        // add the scattering levels, if requested, even if they are all zero
        Array empty;
        if (_recordComponents)
            for (int i = 0; i != _numScatteringLevels; ++i)
            {
                sedNames.push_back(std::to_string(i + 1) + "-times scattered primary flux");
                sedArrays.push_back(_recordTotalOnly ? &empty : &_sed[PrimaryScatteredLevel + i]);
            }

        // construct header comment line
        string header = "# SED at ";
        header += "inclination " + StringUtils::toString(units->oposangle(_inclination)) + " " + units->uposangle();
        header += ", azimuth " + StringUtils::toString(units->oposangle(_azimuth)) + " " + units->uposangle();
        if (_recordPolarization)
        {
            header += ", roll " + StringUtils::toString(units->oposangle(_roll)) + " " + units->uposangle();
        }
        if (_redshift)
        {
            header += ", redshift " + StringUtils::toString(_redshift);
            header += ", luminosity distance " + StringUtils::toString(units->odistance(_luminosityDistance)) + " "
                      + units->udistance();
        }
        else
        {
            header +=
                ", distance " + StringUtils::toString(units->odistance(_luminosityDistance)) + " " + units->udistance();
        }

        // open the file and add the column headers
        TextOutFile sedFile(_parentItem, _instrumentName + "_sed", "SED");
        sedFile.writeLine(header);
        sedFile.addColumn("wavelength; " + units->swavelength(), units->uwavelength());
        for (const string& name : sedNames)
        {
            sedFile.addColumn(name + "; " + units->sfluxdensity(), units->ufluxdensity());
        }

        // write the column data
        for (int ell : Indices(numWavelengths, units->rwavelength()))
        {
            vector<double> values({units->owavelength(_lambdagrid->wavelength(ell))});
            for (const Array* array : sedArrays) values.push_back(array->size() ? (*array)[ell] : 0.);
            sedFile.writeRow(values);
        }
        sedFile.close();

        // output statistics to a seperate file
        if (_recordStatistics)
        {
            // open the file and add the column headers
            TextOutFile statFile(_parentItem, _instrumentName + "_sedstats", "SED statistics");
            statFile.addColumn("wavelength; " + units->swavelength(), units->uwavelength());
            for (int k = 0; k <= maxContributionPower; ++k)
            {
                statFile.addColumn("Sum[w_i**" + std::to_string(k) + "]");
            }
            statFile.writeLine("# --> w_i is luminosity contribution (in W) from i_th launched photon");

            // write the column data
            for (int ell : Indices(numWavelengths, units->rwavelength()))
            {
                vector<double> values({units->owavelength(_lambdagrid->wavelength(ell))});
                for (int k = 0; k <= maxContributionPower; ++k) values.push_back(_wsed[k][ell]);
                statFile.writeRow(values);
            }
            statFile.close();
        }
    }

    // write IFUs to FITS files (one file per IFU)
    if (_includeSurfaceBrightness)
    {
        // Build a list of file names and corresponding pointers to ifu arrays (which may be empty)
        vector<string> ifuNames;
        vector<Array*> ifuArrays;

        // add the total flux; if we didn't record it directly, calculate it now
        ifuNames.push_back("total");
        Array ifuTotal;
        if (_recordTotalOnly)
            ifuArrays.push_back(&_ifu[Total]);
        else
        {
            ifuTotal = _ifu[PrimaryDirect] + _ifu[PrimaryScattered];
            if (_hasMediumEmission) ifuTotal += _ifu[SecondaryDirect] + _ifu[SecondaryScattered];
            ifuArrays.push_back(&ifuTotal);
        }

        // add the flux components, if requested
        if (_recordComponents)
        {
            // add the transparent flux only if it may differ from the total flux
            if (!_recordTotalOnly)
            {
                ifuNames.push_back("transparent");
                ifuArrays.push_back(&_ifu[Transparent]);
            }
            // add the actual components of the total flux (empty arrays will be ignored later on)
            ifuNames.insert(ifuNames.end(), {"primarydirect", "primaryscattered", "secondarytransparent",
                                             "secondarydirect", "secondaryscattered"});
            ifuArrays.insert(ifuArrays.end(),
                             {&_ifu[PrimaryDirect], &_ifu[PrimaryScattered], &_ifu[SecondaryTransparent],
                              &_ifu[SecondaryDirect], &_ifu[SecondaryScattered]});
        }

        // add the polarization components, if requested
        if (_recordPolarization)
        {
            ifuNames.insert(ifuNames.end(), {"stokesQ", "stokesU", "stokesV"});
            ifuArrays.insert(ifuArrays.end(), {&_ifu[TotalQ], &_ifu[TotalU], &_ifu[TotalV]});
        }

        // add the scattering levels, if requested
        if (!_recordTotalOnly)
            for (int i = 0; i != _numScatteringLevels; ++i)
            {
                ifuNames.push_back("primaryscattered" + std::to_string(i + 1));
                ifuArrays.push_back(&_ifu[PrimaryScatteredLevel + i]);
            }

        // copy the wavelength grid in output units
        Array wavegrid(numWavelengths);
        for (int ell = 0; ell != numWavelengths; ++ell)
            wavegrid[ell] = units->owavelength(_lambdagrid->wavelength(ell));

        // reverse the ordering of the wavelength grid and frames if necessary
        if (units->rwavelength())
        {
            NR::reverse(wavegrid);

            // flux frames
            for (auto array : ifuArrays)
                if (array->size()) NR::reverse(*array, _numPixelsInFrame);

            // statistics frames
            for (auto& array : _wifu)
                if (array.size()) NR::reverse(array, _numPixelsInFrame);
        }

        // determine spatial axes values and units
        double incx, incy, cx, cy;
        string unitsxy;
        if (_local)
        {
            // for local instruments, use the metadata provided by the instrument
            if (_quantityXY.empty())
            {
                incx = _incrementX;
                incy = _incrementY;
                cx = _centerX;
                cy = _centerY;
                unitsxy = "1";
            }
            else
            {
                incx = units->out(_quantityXY, _incrementX);
                incy = units->out(_quantityXY, _incrementY);
                cx = units->out(_quantityXY, _centerX);
                cy = units->out(_quantityXY, _centerY);
                unitsxy = units->unit(_quantityXY);
            }
        }
        else
        {
            // for distant instruments, convert to angular sizes
            incx = units->oangle(2. * atan(0.5 * _pixelSizeX / _angularDiameterDistance));
            incy = units->oangle(2. * atan(0.5 * _pixelSizeY / _angularDiameterDistance));
            cx = units->oangle(2. * atan(0.5 * _centerX / _angularDiameterDistance));
            cy = units->oangle(2. * atan(0.5 * _centerY / _angularDiameterDistance));
            unitsxy = units->uangle();
        }

        // determine observer info for distant instruments
        std::unique_ptr<FITSInOut::ObserverInfo> obsInfo;
        if (!_local)
        {
            obsInfo = std::make_unique<FITSInOut::ObserverInfo>();
            obsInfo->inclination = _inclination * (180. / M_PI);
            obsInfo->azimuth = _azimuth * (180. / M_PI);
            obsInfo->roll = _roll * (180. / M_PI);
            obsInfo->redshift = _redshift;
            obsInfo->luminosityDistance = units->odistance(_luminosityDistance);
            obsInfo->angularDiameterDistance = units->odistance(_angularDiameterDistance);
            obsInfo->distanceUnits = units->udistance();
        }

        // output the files (ignoring empty arrays)
        int numFiles = ifuNames.size();
        for (int q = 0; q != numFiles; ++q)
            if (ifuArrays[q]->size())
            {
                string filename = _instrumentName + "_" + ifuNames[q];
                string description = ifuNames[q] + " flux";
                FITSInOut::write(_parentItem, description, filename, *(ifuArrays[q]), units->usurfacebrightness(),
                                 _numPixelsX, _numPixelsY, incx, incy, cx, cy, unitsxy, wavegrid, units->uwavelength(),
                                 obsInfo.get());
            }

        // output statistics to additional files
        if (_recordStatistics)
        {
            // the output files have single-precision floating point numbers with range of only about 10^+-38
            // --> scale the values to a range that has a maximum of 10^+-38 to minimize the number of underflows
            const double WMAX = 1e38;
            Array cs(maxContributionPower);
            for (int k = 1; k <= maxContributionPower; ++k)
            {
                cs[k - 1] = pow(WMAX / _wifu[k].max(), 1. / k);  // inverse of WMAX == c**k w[k].max()
            }
            double c = cs.min();
            double cn = 1.;
            for (int k = 0; k <= maxContributionPower; ++k)
            {
                string filename = _instrumentName + "_stats" + std::to_string(k);
                string description = "sum of contributions to the power of " + std::to_string(k);
                _wifu[k] *= cn;
                FITSInOut::write(_parentItem, description, filename, _wifu[k], "", _numPixelsX, _numPixelsY, incx, incy,
                                 cx, cy, unitsxy, wavegrid, units->uwavelength());
                cn *= c;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

void RayRecorder::recordContributions(ContributionList* contributionList)
{
    // sort the contributions on wavelength and pixel index so that contributions to the same bin are consecutive
    contributionList->sort();
    const vector<Contribution>& contributions = contributionList->contributions();
    size_t numContributions = contributions.size();

    // for SEDs, group contributions on ell index (wavelength bin)
    if (_includeFluxDensity)
    {
        double w = 0;
        for (size_t i = 0; i != numContributions; ++i)
        {
            w += contributions[i].w();
            if (i + 1 == numContributions || contributions[i].ell() != contributions[i + 1].ell())
            {
                int ell = contributions[i].ell();
                double wn = 1.;
                for (int k = 0; k <= maxContributionPower; ++k)
                {
                    LockFree::add(_wsed[k][ell], wn);
                    wn *= w;
                }
                w = 0;
            }
        }
    }

    // for IFUs, group contributions on lell index (wavelength and pixel bins)
    if (_includeSurfaceBrightness)
    {
        double w = 0;
        for (size_t i = 0; i != numContributions; ++i)
        {
            w += contributions[i].w();
            if (i + 1 == numContributions || contributions[i].ell() != contributions[i + 1].ell()
                || contributions[i].l() != contributions[i + 1].l())
            {
                if (contributions[i].l() >= 0)
                {
                    size_t lell = contributions[i].l() + contributions[i].ell() * _numPixelsInFrame;
                    double wn = 1.;
                    for (int k = 0; k <= maxContributionPower; ++k)
                    {
                        LockFree::add(_wifu[k][lell], wn);
                        wn *= w;
                    }
                }
                w = 0;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////
