cmake_minimum_required(VERSION 3.6)
project(MOEA_D)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main_moea.cpp)
add_executable(MOEA_D ${SOURCE_FILES})