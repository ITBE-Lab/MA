# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
# Tips: On Unix with make to get more information: make VERBOSE=1
#
cmake_minimum_required (VERSION 3.8)

# Include contributed (external) projects.


# Define a header only library as taeget
add_library(contrib_CPP_header INTERFACE)
# Include location of header files
target_include_directories(contrib_CPP_header INTERFACE ./json ./intervalTree)

if(WITH_PYTHON)
    if( PythonLibs_FOUND )
        message( "Include pybind11 ..." )
        add_subdirectory (pybind11)
        target_include_directories(contrib_CPP_header SYSTEM INTERFACE ./pybind11_json/include)
    endif()
endif()

# Add pretty printing for tables
# https://github.com/seleznevae/libfort
add_subdirectory(libfort/lib)
