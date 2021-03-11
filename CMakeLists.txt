cmake_minimum_required(VERSION 2.8)

project(cpplex)

set(CMAKE_CXX_FLAGS "-std=c++17 -O3 -Wall -Wextra")

include_directories(${simplex_SOURCE_DIR}/include)

add_executable(solver solver/solver.cc)
add_executable(conicalhull solver/conicalhull.cc)

