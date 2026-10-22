# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library for approximate nearest-neighbour search.
# The *.h header files are taken unmodified from the header-only hnswlib library
# (Hierarchical Navigable Small World graphs) by Yury Malkov et al.
# at https://github.com/nmslib/hnswlib (Apache License 2.0).
#---------------------------------------------------------------------

# set the target name
set(TARGET hnsw)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp" "*.h")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# add SMILE library dependencies
target_link_libraries(${TARGET} fundamentals)
include_directories(../../SMILE/fundamentals)

# adjust C++ compiler flags to our needs
set(NO_EXTRA_WARNINGS true)  # to avoid warnings in the hnswlib code
include("../../SMILE/build/CompilerFlags.cmake")
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(${TARGET} PRIVATE -Wno-type-limits)  # unsigned < 0 comparisons in hnswalg.h
endif()
