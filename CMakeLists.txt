cmake_minimum_required(VERSION 2.8)

project(cpplex)

set(CMAKE_CXX_FLAGS "-std=c++17 -O3")

add_subdirectory(simplex)
include_directories(${simplex_SOURCE_DIR}/include)

add_executable(main solver/solver.cc)
target_link_libraries(main simplex)

