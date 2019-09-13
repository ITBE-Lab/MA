# CMakeList.txt : CMake for minimizerscpp
# Creates the static minimizerscpp library

cmake_minimum_required (VERSION 3.8)

project("minimizerscpp" VERSION 0.1 DESCRIPTION "Minimizer implementation adapted from minimap 2")

# Define the target (library minimizerscpp)
add_library(minimizerscpp STATIC
	src/index.c
	src/bseq.c
	src/kthread.c
	src/misc.c
	src/sketch.c
)

# cmake will add the -fPIC for compilation with GCC
set_property( TARGET minimizerscpp PROPERTY POSITION_INDEPENDENT_CODE ON)

# location of all includes with the kswcpp library
target_include_directories( minimizerscpp PUBLIC inc )

