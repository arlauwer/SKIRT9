/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "IonicCrossSections.hpp"

////////////////////////////////////////////////////////////////////

vector<IonicCrossSections::PhotoAbsorbResource>
IonicCrossSections::loadPhotoAbsorption(SimulationItem* item, const std::unordered_set<int>& ionSet)
{
    return loadPresent<PhotoAbsorbResource, 10>(item, "Ionic_PA.txt", "photo-absorption data", ionSet);
}

////////////////////////////////////////////////////////////////////

vector<IonicCrossSections::FluorescenceResource>
IonicCrossSections::loadFluorescence(SimulationItem* item, const std::unordered_set<int>& ionSet)
{
    return loadPresent<FluorescenceResource, 7>(item, "Ionic_FL.txt", "fluorescence data", ionSet);
}

////////////////////////////////////////////////////////////////////
