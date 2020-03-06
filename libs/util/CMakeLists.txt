# CMakeList.txt : CMake for util
# Creates the static util library

cmake_minimum_required (VERSION 3.8)

project("Utility Library" VERSION 0.1 DESCRIPTION "")

file(GLOB_RECURSE UTIL_SRC
    "src/*.cpp"
) # file


# Define the target (library util)
add_library( util SHARED ${UTIL_SRC} )

# cmake will add the -fPIC for compilation with GCC
set_property( TARGET util PROPERTY POSITION_INDEPENDENT_CODE ON)

# Activate __declspec(dllexport) for DLL export 
# PRIVATE => depending projects will not inherit the EXPORT definition
target_compile_definitions( util PRIVATE DLL_PORT_util=dllexport )

# put the output .dll into the build folder
set_target_properties( util PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )
set_target_properties( util PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

# location of all includes with the util library
target_include_directories( util PUBLIC inc )