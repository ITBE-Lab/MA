# CMakeList.txt : CMake for util
# Creates the static util library

cmake_minimum_required (VERSION 3.8)

project("Utility Library" VERSION 0.1 DESCRIPTION "")

file(GLOB_RECURSE UTIL_SRC
    "src/*.cpp"
) # file


# Define the target (library util)
add_library( util STATIC ${UTIL_SRC} )

# cmake will add the -fPIC for compilation with GCC
set_property( TARGET util PROPERTY POSITION_INDEPENDENT_CODE ON)

# location of all includes with the util library
target_include_directories( util PUBLIC inc )

target_link_libraries ( util LINK_PUBLIC CppSQLite3 )