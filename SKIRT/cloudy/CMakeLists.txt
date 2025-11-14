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
target_link_libraries(${TARGET} fundamentals hnsw utils)
include_directories(../../SMILE/fundamentals ../hnsw ../utils)

# add SKIRT library dependencies
# target_link_libraries(${TARGET} fitsio voro mpi utils tetgen)
# include_directories(../fitsio ../voro ../mpi ../utils)
# include_directories(SYSTEM ../tetgen)  # suppress warnings in tetgen header

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")
