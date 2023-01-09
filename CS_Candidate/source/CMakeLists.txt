# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

set(CMAKE_BUILD_TYPE Release)
#set(BASES_INDICES_LIB ../../BasesToIndices)

# Name of project and that it is C++ only.
project (cs-candidate CXX)


# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for zlib as a dependency for SeqAn.
#find_package (ZLIB)
#find_package (BZip2)
#find_package (OpenMP)


# Load the SeqAn module and fail if not found.
#find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
#include_directories (${SEQAN_INCLUDE_DIRS})
#include_directories(${BASES_INDICES_LIB}/source)

#link_directories(${BASES_INDICES_LIB}/build)

# Add definitions set by find_package (SeqAn).
#add_definitions (${SEQAN_DEFINITIONS})


# Add CXX flags found by find_package (SeqAn).
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")


# Add executable and link against SeqAn dependencies.
add_library(cs-candidate STATIC CS_Candidate.cpp)
#target_link_libraries (cs-candidate ${SEQAN_LIBRARIES})

