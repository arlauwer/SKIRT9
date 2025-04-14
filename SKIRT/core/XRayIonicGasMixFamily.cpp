/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "XRayIonicGasMixFamily.hpp"
#include "Constants.hpp"
#include "StringUtils.hpp"

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

const MaterialMix* XRayIonicGasMixFamily::mix(const Array& parameters)
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

    // use params to make new XRayIonicGasMix
    vector<double> abundances;
    double temperature;

    return new XRayIonicGasMix(this, _ions, boundElectrons, abundances, temperature);
}

////////////////////////////////////////////////////////////////////
