cmake_minimum_required(VERSION 2.8)

project(cpplex)

set(CMAKE_CXX_FLAGS "-std=c++11 -O3")

add_subdirectory(pilal)
add_subdirectory(simplex)

include_directories(${pilal_SOURCE_DIR}/include)
include_directories(${simplex_SOURCE_DIR}/include)

add_executable(solver solver/solver.cc)
target_link_libraries(solver pilal simplex) 

