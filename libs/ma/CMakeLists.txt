cmake_minimum_required (VERSION 3.8)

project("MA Lib" VERSION 1.1.1 DESCRIPTION "Core library of MA")

# Look for boost
# Info: https://feralchicken.wordpress.com/2013/12/07/boost-python-hello-world-example-using-cmake/
# If your boost is not within the standard path, then set BOOST_ROOT appropriately
# set(BOOST_ROOT "/opt/dev/")
# find_package( Boost 1.60.0 COMPONENTS python3 )

# Get all .cpp files in the src folder
file(GLOB_RECURSE MA_SRC
    "src/*.cpp"
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
    add_library( MA SHARED ${MA_SRC} )
endif()

# Activate __declspec(dllexport) for DLL export 
# PRIVATE => depending projects will not inherit the EXPORT definition
target_compile_definitions( MA PRIVATE EXPORT )

# Be verbose with sql queries (in debug mode)
#target_compile_definitions( MA PUBLIC SQL_VERBOSE )

# CMake will add the -fPIC for compilation with GCC
set_property(TARGET MA PROPERTY POSITION_INDEPENDENT_CODE ON)

# Set includes and library dependencies 
target_include_directories(MA PUBLIC inc)
target_link_libraries ( MA LINK_PUBLIC kswcpp CppSQLite3 contrib_CPP_header util )
if( ZLIB_FOUND )
    target_link_libraries ( MA LINK_PUBLIC minimizerscpp )
endif()

# and for each of your dependent executable/library targets:
if( ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0) )
    # std::filesystem is still experimental and for linkage we need separated library
    message( "-- libMA: Filesystem support for GCC below version 8" )
    target_link_libraries( MA PUBLIC stdc++fs )
endif()

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
  ERROR_QUIET
)

# Check if there are any modified files
execute_process(
  COMMAND git status -s
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_STATUS
  RESULT_VARIABLE RESULT_STATUS
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)

if( "${GIT_STATUS}" STREQUAL "" )
    set(GIT_STATUS "")
else()
    set(GIT_STATUS "-D")
endif()

if(NOT RESULT_COMMIT STREQUAL "0" OR NOT RESULT_STATUS STREQUAL "0")
    message( "-- libMA: WARNING cannot infer verison from git (maybe git.exe is not in your environment path?)" )
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

# Integration of test-code, if requested
if(BUILD_TESTING)
    # add tests for all cpp files in the tests folder
    # Get all .cpp files in the test folder
    file(GLOB_RECURSE MA_CPP_TEST_SRC
        "tests/*.cpp"
    ) # file
    foreach(file ${MA_CPP_TEST_SRC})
        get_filename_component(test_name ${file} NAME)
        message( "-- libMA: generating test " ${test_name} )
        add_executable( ${test_name} ${file} )
        target_link_libraries( ${test_name} PUBLIC MA Threads::Threads )
        add_test( ${test_name} ${test_name} )
    endforeach()

    # add tests for all python files in the tests folder
    if( PythonLibs_FOUND )
        find_package(PythonInterp ${requested_python_version})
        if(PythonInterp_FOUND)
            # Get all .py files in the test folder
            file(GLOB_RECURSE MA_PY_TEST_SRC
                "tests/*.py"
            ) # file
            foreach(file ${MA_PY_TEST_SRC})
                get_filename_component(test_name ${file} NAME)
                message( "-- libMA: generating test " ${test_name} )
                add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} ${file} )
            endforeach()
        else()
            message("-- libMA: WARNING failed to generate python tests because no interpreter was found")
        endif()
    endif()
endif()