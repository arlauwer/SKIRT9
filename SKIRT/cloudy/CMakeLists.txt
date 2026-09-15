# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library for the SKIRT Cloudy wrapper
# ------------------------------------------------------------------

# set the target name
set(TARGET cloudy)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# add SMILE library dependencies
target_link_libraries(${TARGET} schema fundamentals)
include_directories(../../SMILE/schema ../../SMILE/fundamentals)

# add required dependencies
target_link_libraries(${TARGET} fundamentals hnsw utils tools wavelength)
include_directories(../../SMILE/fundamentals ../hnsw ../utils ../tools ../wavelength)

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")
