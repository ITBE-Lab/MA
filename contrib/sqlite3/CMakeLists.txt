# CMakeList.txt : CMake for sqlite3
# Creates the static sqlite3 library
cmake_minimum_required(VERSION 3.8)
project("sqlite3" VERSION 0.1 DESCRIPTION "C sqlite library")

add_definitions(-DSQLITE_ENABLE_RTREE)
add_definitions(-DSQLITE_ENABLE_FTS4)
add_definitions(-DSQLITE_ENABLE_FTS5)
add_definitions(-DSQLITE_ENABLE_JSON1)
add_definitions(-DSQLITE_ENABLE_RBU)
add_definitions(-DSQLITE_ENABLE_STAT4)

# Uncomment this for single-threaded variant (faster)
#add_definitions(-DSQLITE_THREADSAFE=0)


# Define the target (library sqlite3)
add_library(sqlite3 STATIC
	src/sqlite3.c
)

# cmake will add the -fPIC for compilation with GCC
set_property( TARGET sqlite3 PROPERTY POSITION_INDEPENDENT_CODE ON)

# location of all includes with the sqlite3 library
target_include_directories( sqlite3 PUBLIC inc )
