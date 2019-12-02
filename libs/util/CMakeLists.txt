# CMakeList.txt : CMake for util
# Creates the static util library

cmake_minimum_required (VERSION 3.8)

project("Utility Library" VERSION 0.1 DESCRIPTION "")

# Define a header only library as taeget
add_library(util INTERFACE)
# Include location of header files
target_include_directories(util INTERFACE inc)