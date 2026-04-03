# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a wrapper for the (header only) hnswlib library
# ------------------------------------------------------------------

# set the target name
set(TARGET hnsw)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# add required dependencies
target_link_libraries(${TARGET} fundamentals utils)
include_directories(../../SMILE/fundamentals ../utils)

# adjust C++ compiler flags to our needs
set(NO_EXTRA_WARNINGS true)  # to avoid warnings in the hnswlib code
include("../../SMILE/build/CompilerFlags.cmake")