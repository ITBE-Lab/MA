# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
# Tips: On Unix with make to get more information: make VERBOSE=1
#
cmake_minimum_required (VERSION 3.8)

# Include sub-projects.
message( "Include sqlite3 ..." )
add_subdirectory (sqlite3)

message( "Include CppSQLite3 ..." )
add_subdirectory (CppSQLite3)
