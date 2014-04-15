cmake_minimum_required(VERSION 2.8)

project(pilal)

set(CMAKE_CXX_FLAGS "-std=c++11 -O3")

include_directories(${pilal_SOURCE_DIR}/include)
add_library(pilal src/matrix.cc)

