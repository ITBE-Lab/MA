# CMakeList.txt : CMake for CppSQLite3
# Creates the static CppSQLite3 library

cmake_minimum_required (VERSION 3.8)

project("CppSQLite3" VERSION 0.1 DESCRIPTION "Cpp Sqlite Library")

# Define the target (library CppSQLite3)
add_library(CppSQLite3 STATIC
	src/CppSQLite3.cpp
)

# cmake will add the -fPIC for compilation with GCC
set_property( TARGET CppSQLite3 PROPERTY POSITION_INDEPENDENT_CODE ON)


target_link_libraries ( CppSQLite3 LINK_PUBLIC sqlite3 )

# location of all includes with the CppSQLite3 library
target_include_directories( CppSQLite3 PUBLIC inc )
