# CMakeList.txt : Top-level CMake project file for the Modular Aligner MA
# Tips: 
#	- For getting more information with make: make VERBOSE=1
#	- Call cmake with "-DCMAKE_BUILD_TYPE=Debug" for switching to the debug mode

cmake_minimum_required (VERSION 3.8)

project("MA" VERSION 0.1 DESCRIPTION "MA - The Modular Aligner")

# optional features
option(WITH_PYTHON "Build python integration of libMA." ON)
option(WITH_ZLIB "Build gz file readers using zlib." ON)
option(WITH_DOCS "Build doxygen documentation." ON)
option(WITH_GUI "Build maGUI" ON)

# assert that the compiler knows c++17
if ( NOT "cxx_std_17" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
   message( FATAL_ERROR "Your compiler does not support C++17, which is required for compiling MA. "
                        "(Updating your compiler will likely fix this issue.)" )
endif()

# Collect all arguments of the CMake call (just for info purposes ...)
get_cmake_property(CACHE_VARS CACHE_VARIABLES)
foreach(CACHE_VAR ${CACHE_VARS})
  get_property(CACHE_VAR_HELPSTRING CACHE ${CACHE_VAR} PROPERTY HELPSTRING)
  if(CACHE_VAR_HELPSTRING STREQUAL "No help, variable specified on the command line.")
    get_property(CACHE_VAR_TYPE CACHE ${CACHE_VAR} PROPERTY TYPE)
    if(CACHE_VAR_TYPE STREQUAL "UNINITIALIZED")
      set(CACHE_VAR_TYPE)
    else()
      set(CACHE_VAR_TYPE :${CACHE_VAR_TYPE})
    endif()
    set(CMAKE_ARGS "${CMAKE_ARGS} -D${CACHE_VAR}${CACHE_VAR_TYPE}=\"${${CACHE_VAR}}\"")
  endif()
endforeach()

# Print info about the current call/environment of CMake
# message("// CMAKE_ARGS: ${CMAKE_ARGS}")
message("// CMAKE_BUILD_TYPE = \"${CMAKE_BUILD_TYPE}\"")
# message("// CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
# message("// CMAKE_CXX_FLAGS_DEBUG = ${CMAKE_CXX_FLAGS_DEBUG}")
# message("// CMAKE_CXX_FLAGS_RELEASE = ${CMAKE_CXX_FLAGS_RELEASE}")
# message("// CMAKE_CXX_FLAGS_RELWITHDEBINFO = ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
# message("// CMAKE_CXX_FLAGS_MINSIZEREL = ${CMAKE_CXX_FLAGS_MINSIZEREL}")

if( CMAKE_BUILD_TYPE STREQUAL "" )
	# Switch to the release mode, so that we get full optimization
	message( "Build type set to \"Release\"" )
	set( CMAKE_BUILD_TYPE "Release" )
endif()

# Currently, the debug_level settings effect all projects (incl. kswcpp)
if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
	# Set the debug level 
	message( "Set debug level to 1" )
	add_definitions( "-DDEBUG_LEVEL=1" )
	if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
		# In the debug mode the code comprises unsed variables.
		# Corresponding warning will be switched off.
		message( "Add -Wno-unused-variable for gcc")
		add_compile_options( "-Wno-unused-variable" )
	elseif( MSVC )
		# FOr MSVC with have to adapt the RelWithDebInfo setting in order to avoid linking problems
		set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "/MDd /Zi /O2 /Ob2")
		message("CMAKE_CXX_FLAGS_RELWITHDEBINFO changed to ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
	endif()
endif()

# Testing with cmake
# https://stackoverflow.com/questions/50468620/what-does-enable-testing-do-in-cmake
option(ENABLE_TESTS "Enable tests" ON)
if (${ENABLE_TESTS})
    include(CTest)
endif()

# Add additional general compile flags
# Note: Some flags are already preset by the generator depending on the build-type
# TODO: use set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -foo" )
if( NOT MSVC )
	# Note: Adding "-Wextra" delivers more warnings
    # "-Werror" "-Wall" : All warnings enabled and warning are treated as errors
    add_compile_options( "-Wall" )
    if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
        #add_compile_options( "-Werror" )
    endif()
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        add_compile_options( "-stdlib=libc++" )
    endif()
endif()

if(WITH_PYTHON)
    # Look for Python
    set(requested_python_version 3.5)
    find_package( PythonLibs ${requested_python_version} )
    if( PythonLibs_FOUND )
        message( "Found Python Libs")
    endif()
endif()

if(WITH_ZLIB)
    # Look for zlib
    find_package( ZLIB )
    if( ZLIB_FOUND )
        message( "Found zlib" )
    endif()
endif()

if(WITH_DOCS)
    # check if Doxygen is installed
    find_package(Doxygen)
    if (DOXYGEN_FOUND)

        set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doxygen.config)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config.out)
        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
        # note the option ALL which allows to build the docs together with the application
        add_custom_target( docs
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM )
    else (DOXYGEN_FOUND)
        message("Doxygen needs to be installed in order to generate the documentation")
    endif (DOXYGEN_FOUND)
endif()

# Required for linking via libpthreads on UNIX systems
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# Use Cpp 17 standard
if (NOT MSVC)
    set(PYBIND11_CPP_STANDARD -std=c++17)
else()
    set(PYBIND11_CPP_STANDARD /std:c++17)
endif()
set(CMAKE_CXX_STANDARD 17)
# Build all required contributor libraries
add_subdirectory ( contrib )
# Build all required libraries
add_subdirectory ( libs )
# Build wxWidgets GUI
if(WITH_GUI)
    add_subdirectory ( gui )
endif()

# Build the MA executable itself
add_executable ( maCMD cmdMa.cpp )
# Link cmdMA
target_link_libraries ( maCMD PUBLIC MA Threads::Threads )
