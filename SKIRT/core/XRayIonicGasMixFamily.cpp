/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "XRayIonicGasMixFamily.hpp"
#include "Constants.hpp"
#include "StringUtils.hpp"

XRayIonicGasMixFamily::~XRayIonicGasMixFamily()
{
    for (XRayIonicGasMix* mix : _mixes)
    {
        delete mix;
    }
}

void XRayIonicGasMixFamily::setupSelfBefore()
{
    MaterialMixFamily::setupSelfBefore();

    // read ions
    for (string ion : StringUtils::split(ions(), ","))
    {
        ion = StringUtils::squeeze(ion);
        _ionNames.push_back(ion);
    }
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> XRayIonicGasMixFamily::parameterInfo() const
{
    vector<SnapshotParameter> descriptors;

    for (string ionName : _ionNames)
    {
        descriptors.push_back(SnapshotParameter::custom(ionName));
    }

    return descriptors;
}

////////////////////////////////////////////////////////////////////

const MaterialMix* XRayIonicGasMixFamily::mix(double /*Z*/, double T, const Array& parameters)
{
    XRayIonicGasMix::BoundElectrons boundElectrons;
    switch (_scatterBoundElectrons)
    {
        case BoundElectrons::None: boundElectrons = XRayIonicGasMix::BoundElectrons::None; break;
        case BoundElectrons::Free: boundElectrons = XRayIonicGasMix::BoundElectrons::Free; break;
        case BoundElectrons::FreeWithPolarization:
            boundElectrons = XRayIonicGasMix::BoundElectrons::FreeWithPolarization;
            break;
        case BoundElectrons::Good: boundElectrons = XRayIonicGasMix::BoundElectrons::Good; break;
        case BoundElectrons::Exact: boundElectrons = XRayIonicGasMix::BoundElectrons::Exact; break;
    }

    // convert Array to vector
    vector<double> abundances(begin(parameters), end(parameters));

    // look for duplicates
    for (const XRayIonicGasMix* mix : _mixes)
    {
        if (mix->abundances() == abundances) return mix;
    }

    XRayIonicGasMix* mix = new XRayIonicGasMix(this, _ions, boundElectrons, abundances, T, true);
    _mixes.push_back(mix);

    return mix;
}

////////////////////////////////////////////////////////////////////

const MaterialMix* XRayIonicGasMixFamily::mix()
{
    if (_mixes.size() > 0) return _mixes[0];

    // create a default mix for the Configuration (might remove later)
    if (!_defaultMix)
    {
        _defaultMix = new XRayIonicGasMix(this, _ions, XRayIonicGasMix::BoundElectrons::None,
                                          vector<double>(_ionNames.size(), 0.), 0., false);
    }
    return _defaultMix;
}
