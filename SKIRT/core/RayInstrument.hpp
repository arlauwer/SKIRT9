/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RAYINSTRUMENT_HPP
#define RAYINSTRUMENT_HPP

#include "Direction.hpp"
#include "Position.hpp"
#include "SimulationItem.hpp"
#include "WavelengthGrid.hpp"
#include "MediumSystem.hpp"
class FluxRecorder;
class PhotonPacket;

////////////////////////////////////////////////////////////////////

/** RayInstrument is an abstract class representing instruments to collect the photon packets during a
    Monte Carlo simulation. Subclasses implement instruments that can vary in type of projection
    (e.g. parallel or perspective) and in what is recorded (e.g. SED or full integral field data
    cube). This top-level abstract class offers a generic interface for receiving photon packets
    from the simulation. It also includes facilities for configuring user properties that are
    common to all instruments, such as which flux contributions need to be recorded. A wavelength
    grid is established either by specifying a grid for this instrument specifically, or by
    defaulting to the common grid specified for the instrument system. */
class RayInstrument : public SimulationItem
{
    ITEM_ABSTRACT(RayInstrument, SimulationItem, "an instrument")

        PROPERTY_STRING(instrumentName, "the name for this instrument")

        PROPERTY_ITEM(wavelengthGrid, WavelengthGrid, "the wavelength grid for this instrument")
        ATTRIBUTE_RELEVANT_IF(wavelengthGrid, "Panchromatic")
        ATTRIBUTE_REQUIRED_IF(wavelengthGrid, "!DefaultRayInstrumentWavelengthGrid")
        ATTRIBUTE_DISPLAYED_IF(wavelengthGrid, "Level2")

        ATTRIBUTE_SUB_PROPERTIES_HERE(RayInstrument)

        PROPERTY_BOOL(recordComponents, "record flux components separately")
        ATTRIBUTE_DEFAULT_VALUE(recordComponents, "false")

        PROPERTY_INT(numScatteringLevels, "the number of individually recorded scattering levels")
        ATTRIBUTE_MIN_VALUE(numScatteringLevels, "0")
        ATTRIBUTE_MAX_VALUE(numScatteringLevels, "99")
        ATTRIBUTE_DEFAULT_VALUE(numScatteringLevels, "0")
        ATTRIBUTE_RELEVANT_IF(numScatteringLevels, "recordComponents")
        ATTRIBUTE_DISPLAYED_IF(numScatteringLevels, "Level2")

        PROPERTY_BOOL(recordPolarization, "record polarization (Stokes vector elements)")
        ATTRIBUTE_DEFAULT_VALUE(recordPolarization, "false")
        ATTRIBUTE_DISPLAYED_IF(recordPolarization, "Level2")

        PROPERTY_BOOL(recordStatistics, "record information for calculating statistical properties")
        ATTRIBUTE_DEFAULT_VALUE(recordStatistics, "false")
        ATTRIBUTE_DISPLAYED_IF(recordStatistics, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function performs setup for the instrument. It establishes the wavelength grid for the
        instrument: if a grid is specified for the instrument, that grid is used. If not, the
        default wavelength grid specified for the instrument system is used instead. If neither of
        these grids are specified, the function throws a fatal error.

        The function also creates and partially configures the FluxRecorder instance for this
        instrument, passing it the values of the user properties offered by this class and some
        extra information on the simulation. The setupSelfBefore() function of each subclass is
        expected to augment the configuration by calling the includeFluxDensity() and/or
        includeSurfaceBrightness() functions. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the instrument name as human-readable name for the simulation item,
        so that it can be used in log messages to identify the instrument and differentiate it from
        other instruments. */
    string itemName() const override;

    /** This function returns the wavelength grid for the instrument as determined during setup,
        i.e. either the grid specified for this instrument or the default grid specified for the
        instrument system. After setup has completed, the function never returns a nulltpr because
        setupSelfBefore() throws a fatal error if neither of these grids are specified. */
    const WavelengthGrid* instrumentWavelengthGrid() const { return _instrumentWavelengthGrid; }

    virtual void rayTrace() = 0;

    /** This function calibrates the instrument and outputs the recorded contents to a set of
        files. It simply calls the corresponding function of the FluxRecorder instance associated
        with this instrument. */
    void write();

    /** This function returns true if the receiving instrument has the same observer type, position
        and viewing direction as the preceding instrument in the instrument system. This
        information is determined and cached by the determineSameObserverAsPreceding() function,
        which is called by the RayInstrumentSystem during setup. */
    bool isSameObserverAsPreceding() const { return _isSameObserverAsPreceding; }

protected:
    /** This function sets the "isSameObserverAsPreceding" flag to true. By default (i.e. if this
        function is never invoked) the flag is set to false. This function is intended for use from
        the determineSameObserverAsPreceding() function implementation in subclasses. */
    void setSameObserverAsPreceding() { _isSameObserverAsPreceding = true; }

private:
    void trace(Position bfr, Direction bfk);

    //=========== Functions to be implemented in subclass ===========

public:
    /** This function determines whether the specified instrument has the same observer type,
        position and viewing direction as the receiving instrument, and if so, calls the
        setSameObserverAsPreceding() function to remember the fact. The function is invoked by the
        RayInstrumentSystem during setup. The implementation must be provided in a subclass. */
    virtual void determineSameObserverAsPreceding(const RayInstrument* precedingRayInstrument) = 0;

    /** This function returns the direction towards the observer, expressed in model coordinates,
        given the photon packet's launching position. The implementation must be provided in a
        subclass. */
    virtual Direction bfkobs(const Position& bfr) const = 0;

    /** This function returns the direction along the positive y-axis of the instrument frame,
        expressed in model coordinates, given the photon packet's launching position. The
        implementation must be provided in a subclass. */
    virtual Direction bfky(const Position& bfr) const = 0;

    //======================== Data Members =======================

private:
    const WavelengthGrid* _instrumentWavelengthGrid{nullptr};
    bool _isSameObserverAsPreceding{false};
    MediumSystem* _ms;
    SpatialGrid* _grid;
};

////////////////////////////////////////////////////////////////////

#endif
