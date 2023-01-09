# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

set(CMAKE_BUILD_TYPE Release)
set(BASES_INDICES_LIB ../../BasesToIndices/)

# Name of project and that it is C++ only.
project (CS_Consensus CXX)


# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

include_directories(${BASES_INDICES_LIB}/source)

link_directories(${BASES_INDICES_LIB}/build)

add_library(CS_Consensus STATIC CS_Consensus.cpp)
target_link_libraries (CS_Consensus BasesIndices)

