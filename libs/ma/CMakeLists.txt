cmake_minimum_required (VERSION 3.8)

project("MA Lib" VERSION 1.0.1 DESCRIPTION "Core library of MA")

# Look for boost
# Info: https://feralchicken.wordpress.com/2013/12/07/boost-python-hello-world-example-using-cmake/
# If your boost is not within the standard path, then set BOOST_ROOT appropriately
# set(BOOST_ROOT "/opt/dev/")
# find_package( Boost 1.60.0 COMPONENTS python3 )

# Get all .cpp files in the src folder
file(GLOB MA_SRC
    "src/*/*.cpp"
) # file

# Define library MA as target
if( PythonLibs_FOUND )
    # Build shared library with Python support and surpress warnings form within Pybind11
    message( "-- libMA: Build shared library with python support via local pybind11." )
    pybind11_add_module( MA SHARED SYSTEM ${MA_SRC} )
    target_compile_definitions( MA PRIVATE WITH_PYTHON )
else()
    # Build the plain shared library without Python support
    message( "-- libMA: Build shared library without python support." )
    add_library(MA SHARED ${MA_SRC} ) 
endif()

# Activate __declspec(dllexport) for DLL export 
# PRIVATE => depending projects will not inherit the EXPORT definition
target_compile_definitions( MA PRIVATE EXPORT )

# CMake will add the -fPIC for compilation with GCC
set_property(TARGET MA PROPERTY POSITION_INDEPENDENT_CODE ON)

# Set includes and library dependencies 
target_include_directories(MA PUBLIC inc ../kswcpp/inc)
target_link_libraries ( MA LINK_PUBLIC kswcpp )

# if there is a zlib
if( ZLIB_FOUND )
    # Activate code that relies on zlib by defining WITH_ZLIB
    target_compile_definitions( MA PUBLIC WITH_ZLIB )
    # Add include directories of zlib
    target_include_directories(MA PUBLIC ${ZLIB_INCLUDE_DIRS})
    # link to the zlib libraries
    target_link_libraries ( MA LINK_PUBLIC ${ZLIB_LIBRARIES} )
    message( "-- libMA: Building shared library components that rely on zlib." )
else()
    message( "-- libMA: Not building shared library components that rely on zlib." )
endif()

# Generate a version of MA based on the current commit

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  RESULT_VARIABLE RESULT_COMMIT
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Check if there are any modified files
execute_process(
  COMMAND git status -s
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_STATUS
  RESULT_VARIABLE RESULT_STATUS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if( "${GIT_STATUS}" STREQUAL "" )
    set(GIT_STATUS "")
else()
    set(GIT_STATUS "-D")
endif()

if(NOT RESULT_COMMIT STREQUAL "0" OR NOT RESULT_STATUS STREQUAL "0")
    message( "-- libMA: WARNING cannot inferr verison from git (maybe git.exe is not in your environment path?)" )
    set(MA_VERSION "Unknown")
else()
    set(MA_VERSION "${PROJECT_VERSION}-${GIT_COMMIT_HASH}${GIT_STATUS}")
    message( "-- libMA: Version = ${MA_VERSION}" )
endif()

configure_file( # @todo so far this only triggers during configure not at build time
  inc/util/version.h.in
  ${CMAKE_BINARY_DIR}/generated/version.h
)
target_include_directories(MA PUBLIC ${CMAKE_BINARY_DIR}/generated)

# build the shared library in the folder of the .exe 
set_target_properties(MA PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )
set_target_properties(MA PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

# additionally link to python if available
if( PythonLibs_FOUND )
    target_link_libraries ( MA PUBLIC ${PYTHON_LIBRARIES} )
    # if we are under windows and we want to support python we need to do two things:
    if( MSVC )
        # change the name of the library from MA to libMA (under linux lib is used as prefix automatically)
        set_target_properties( MA PROPERTIES OUTPUT_NAME "libMA" )
        # change the suffix to .pdy so that python recognized the file as module 
        # (this does not effect the executable; it can still load libMA.pyd as dynamic linked library)
        set_target_properties( MA PROPERTIES SUFFIX ".pyd" )
    endif()
endif()