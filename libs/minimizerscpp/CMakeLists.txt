# CMakeList.txt : CMake for minimizerscpp
# Creates the static minimizerscpp library

cmake_minimum_required (VERSION 3.8)

project("minimizerscpp" VERSION 0.1 DESCRIPTION "Minimizer implementation adapted from minimap 2")

# Define the target (library minimizerscpp)
add_library(minimizerscpp STATIC
	src/index.cpp
	src/bseq.cpp
	src/kthread.cpp
	src/misc.cpp
	src/sketch.cpp
	src/map.cpp
	src/sdust.cpp
	src/kalloc.cpp
)

# if there is a zlib
if( ZLIB_FOUND )
    # Activate code that relies on zlib by defining WITH_ZLIB
    target_compile_definitions( minimizerscpp PUBLIC WITH_ZLIB )
    # Add include directories of zlib
    target_include_directories( minimizerscpp SYSTEM PUBLIC ${ZLIB_INCLUDE_DIRS})
    # link to the zlib libraries
    target_link_libraries ( minimizerscpp LINK_PUBLIC ${ZLIB_LIBRARIES} )
else()
    message( FATAL_ERROR "-- libminimizerscpp: zlib is required." )
endif()

# cmake will add the -fPIC for compilation with GCC
set_property( TARGET minimizerscpp PROPERTY POSITION_INDEPENDENT_CODE ON)

# location of all includes with the kswcpp library
target_include_directories( minimizerscpp PUBLIC inc )

