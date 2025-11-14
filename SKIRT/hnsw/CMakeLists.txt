# set the target name
set(TARGET hnsw)

add_library(hnsw INTERFACE)

# target_include_directories(hnsw INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# adjust C++ compiler flags to our needs
# set(NO_EXTRA_WARNINGS true)
# include("../../SMILE/build/CompilerFlags.cmake")
