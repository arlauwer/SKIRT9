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

////////////////////////////////////////////////////////////////////

void XRayIonicGasMixFamily::setupSelfBefore()
{
    MaterialMixFamily::setupSelfBefore();

    setup();
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
    // convert Array to vector
    vector<double> abundances(begin(parameters), end(parameters));

    // look for duplicates
    for (const XRayIonicGasMix* mix : _mixes)
    {
        if (mix->abundances() == abundances) return mix;
    }

    XRayIonicGasMix* mix = new XRayIonicGasMix(this, _ions, abundances, T, _boundElectrons, true);
    _mixes.push_back(mix);

    return mix;
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMixFamily::setup()
{
    // read ions (temp fix)
    if (_ionNames.size() == 0)
    {
        for (string ion : StringUtils::split(ions(), ","))
        {
            ion = StringUtils::squeeze(ion);
            _ionNames.push_back(ion);
        }
    }

    // convert enum
    switch (scatterBoundElectrons())
    {
        case BoundElectrons::None: _boundElectrons = XRayIonicGasMix::BoundElectrons::None; break;
        case BoundElectrons::Free: _boundElectrons = XRayIonicGasMix::BoundElectrons::Free; break;
        case BoundElectrons::FreeWithPolarization:
            _boundElectrons = XRayIonicGasMix::BoundElectrons::FreeWithPolarization;
            break;
        case BoundElectrons::Good: _boundElectrons = XRayIonicGasMix::BoundElectrons::Good; break;
        case BoundElectrons::Exact: _boundElectrons = XRayIonicGasMix::BoundElectrons::Exact; break;
    }
}

////////////////////////////////////////////////////////////////////

const MaterialMix* XRayIonicGasMixFamily::mix()
{
    if (_mixes.size() > 0) return _mixes[0];

    // create a default mix for the Configuration (might remove later)
    if (!_defaultMix)
    {
        setup();

        vector<double> abundances(_ionNames.size(), 0.);
        _defaultMix = new XRayIonicGasMix(this, _ions, abundances, 0., _boundElectrons, true);
    }
    return _defaultMix;
}
