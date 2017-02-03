cmake_minimum_required(VERSION 3.6)
project(diamond)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

#set(SOURCE_FILES main.cpp)
set(SOURCE_FILES correct.cpp)
add_executable(diamond ${SOURCE_FILES})