# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library that runs the external photoionization code Cloudy
# (https://gitlab.nublado.org/cloudy/cloudy).
# ------------------------------------------------------------------

# set the target name
set(TARGET cloudy)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# add SMILE library dependencies
target_link_libraries(${TARGET} fundamentals)
include_directories(../../SMILE/fundamentals)

# add direct SKIRT library dependencies
target_link_libraries(${TARGET} hnsw)
include_directories(../hnsw)

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")
