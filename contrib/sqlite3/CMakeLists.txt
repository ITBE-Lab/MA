# CMakeList.txt : CMake for sqlite3
# Creates the static sqlite3 library
cmake_minimum_required(VERSION 3.8)
project("sqlite3" VERSION 3.25.3 DESCRIPTION "C sqlite library")

add_definitions( -DSQLITE_ENABLE_RTREE ) # Include R-Trees for range queries
# increase the default page size for insert performance
# https://www.sqlite.org/pragma.html#pragma_page_size
# https://stackoverflow.com/questions/1711631/improve-insert-per-second-performance-of-sqlite
add_definitions( -DSQLITE_DEFAULT_PAGE_SIZE=65536 )

# Uncomment this to allow users to create special tables with a built-in full-text index
# add_definitions( -DSQLITE_ENABLE_FTS4 )
# add_definitions( -DSQLITE_ENABLE_FTS5 )

# Uncomment this for getting additonal JSON support
# add_definitions( -DSQLITE_ENABLE_JSON1 )

# Define the target (library sqlite3)
add_library(sqlite3 STATIC
	src/sqlite3.c
)

if( CMAKE_BUILD_TYPE STREQUAL "Release" )
    message( "Sqlite3 is compiled using SQLITE_THREADSAFE=0" )
    add_definitions( -DSQLITE_THREADSAFE=0 )
endif()

# cmake will add the -fPIC for compilation with GCC
set_property( TARGET sqlite3 PROPERTY POSITION_INDEPENDENT_CODE ON )

# location of all includes with the sqlite3 library
target_include_directories( sqlite3 PUBLIC inc )
