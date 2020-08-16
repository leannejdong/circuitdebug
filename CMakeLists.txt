cmake_minimum_required(VERSION 3.17)
project(circuitdebug)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG")

add_executable(circuitdebug algorithm.cpp)