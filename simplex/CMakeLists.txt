cmake_minimum_required(VERSION 2.8)

project(simplex)

set(CMAKE_CXX_FLAGS "-std=c++11 -O3")

include_directories(${simplex_SOURCE_DIR}/include)
include_directories(${pilal_SOURCE_DIR}/include)


add_library(simplex src/columnset.cc src/constraint.cc src/objectivefunction.cc src/simplex.cc src/variable.cc)

