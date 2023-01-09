# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

set(CMAKE_BUILD_TYPE Release)

# Name of project and that it is C++ only.
project (baseToIndice CXX)


# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.


# Add executable and link against SeqAn dependencies.
add_library(BasesIndices STATIC BasesIndices.cpp)


