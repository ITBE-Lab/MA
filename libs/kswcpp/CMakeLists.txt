# CMakeList.txt : CMake for kswcpp
# Creates the static kswcpp library

cmake_minimum_required (VERSION 3.8)

project("kswcpp" VERSION 0.1 DESCRIPTION "Needleman-Wunsch SSE2 implementation")

# Define the target (library kswcpp)
add_library(kswcpp STATIC
	src/cpu_info.cpp
	src/kswcpp_sse_xx.cpp
	src/kswcpp_avx2.cpp
)

# cmake will add the -fPIC for compilation with GCC
set_property( TARGET kswcpp PROPERTY POSITION_INDEPENDENT_CODE ON)

# File-depend compile-flags:
# https://stackoverflow.com/questions/20244800/appending-compiler-flags-to-a-file-with-cmake
if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
	# using GCC
	set_property( SOURCE "src/kswcpp_sse_xx.cpp" APPEND_STRING PROPERTY COMPILE_FLAGS "-msse4.1" )
	set_property( SOURCE "src/kswcpp_avx2.cpp" APPEND_STRING PROPERTY COMPILE_FLAGS "-mavx2" )
elseif( MSVC )
	# using Visual Studio C++
	add_definitions(-DFOO) # Workaround (see https://github.com/Microsoft/VSLinux/issues/288)
	set_property( SOURCE src/kswcpp_avx2.cpp APPEND_STRING PROPERTY COMPILE_FLAGS "/arch:AVX2" )
else()
	message( FATAL_ERROR "kswcpp does not support your compiler now" )
endif()

# location of all includes with the kswcpp library
target_include_directories( kswcpp PUBLIC inc )

