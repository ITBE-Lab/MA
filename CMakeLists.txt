# CMakeList.txt : Top-level CMake project file for the Modular Aligner MA
# Tips: 
#	- For getting more information with make: make VERBOSE=1
#	- Call cmake with "-DCMAKE_BUILD_TYPE=Debug" for switching to the debug mode

cmake_minimum_required (VERSION 3.8)

project("MA" VERSION 0.1 DESCRIPTION "MA - The Modular Aligner")


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
message("// CMAKE_ARGS: ${CMAKE_ARGS}")
message("// CMAKE_BUILD_TYPE = \"${CMAKE_BUILD_TYPE}\"")
message("// CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
message("// CMAKE_CXX_FLAGS_DEBUG = ${CMAKE_CXX_FLAGS_DEBUG}")
message("// CMAKE_CXX_FLAGS_RELEASE = ${CMAKE_CXX_FLAGS_RELEASE}")
message("// CMAKE_CXX_FLAGS_RELWITHDEBINFO = ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
message("// CMAKE_CXX_FLAGS_MINSIZEREL = ${CMAKE_CXX_FLAGS_MINSIZEREL}")

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
		# In the debuf mode the code comprises unsed variables.
		# Corresponding warning will be switched off.
		message( "Add -Wno-unused-variable for gcc")
		add_compile_options( "-Wno-unused-variable" )
	elseif( MSVC )
		# FOr MSVC with have to adapt the RelWithDebInfo setting in order to avoid linking problems
		set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "/MDd /Zi /O2 /Ob2")
		message("CMAKE_CXX_FLAGS_RELWITHDEBINFO changed to ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
	endif()
endif()

# Add additional general compile flags
# Note: Some flags are already preset by the generator depending on the build-type
# TODO: use set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -foo" )
if( NOT MSVC )
	# Note: Adding "-Wextra" delivers more warnings 
    add_compile_options( "-Werror" "-Wall" ) 
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        add_compile_options( "-stdlib=libc++" )
    endif()
endif()

# Look for Python
find_package( PythonLibs 3.5 )
if( PythonLibs_FOUND )
	message( "Found Python Libs")
endif()

# Look for zlib
find_package( ZLIB )
if( ZLIB_FOUND )
	message( "Found zlib" )
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

# Build the MA excecutable itself
add_executable ( maCMD cmdMa.cpp )
# Link cmdMA
target_link_libraries ( maCMD PUBLIC MA Threads::Threads )
