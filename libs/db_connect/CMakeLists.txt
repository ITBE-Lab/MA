cmake_minimum_required (VERSION 3.8)

project("db_connect" VERSION 0.1 DESCRIPTION "Database connection library")

INCLUDE(cmake/FindMySQL.cmake)

if( MySQL_FOUND )
        message( "Found MySQL includes and libraries")
		# location of all includes with the kswcpp library
		#target_include_directories( db_connect PUBLIC ${MySQL_INCLUDE_DIRS} )
		#target_link_libraries ( db_connect PUBLIC ${MySQL_LIBRARIES} )
endif()

# Define the target (library CppSQLite3)
add_library(dbSQL STATIC
	connectors/CppSQLite3.cpp db_sql.cpp
)
# cmake will add the -fPIC for compilation with GCC
set_property( TARGET dbSQL PROPERTY POSITION_INDEPENDENT_CODE ON)

# Insert the sqlite3 object file
target_link_libraries ( dbSQL LINK_PUBLIC sqlite3 ${MYSQL_LIBRARIES} )

# location of all includes 
target_include_directories( dbSQL PUBLIC connectors . ${MYSQL_INCLUDE_DIRS})

# # Define a header only library as taeget
# add_library(dbSQL INTERFACE)
# # Include location of header files
# target_include_directories(dbSQL INTERFACE ./)

message( "Makefile.txt - Include:" ${MYSQL_INCLUDE_DIRS} " Library: " ${MYSQL_LIBRARIES} )

# Integration of test-code, if requested
if(BUILD_TESTING)
    # add tests for all cpp files in the tests folder
    # Get all .cpp files in the test folder
    file(GLOB MA_CPP_TEST_SRC
        "tests/*.cpp"
    ) # file
    foreach(file ${MA_CPP_TEST_SRC})
        get_filename_component(test_name ${file} NAME)
        message( "-- MA Database Connector: generating test " ${test_name} )
        add_executable( ${test_name} ${file} )
		target_include_directories( ${test_name} PUBLIC ${MYSQL_INCLUDE_DIRS} connectors   )
        target_link_libraries( ${test_name} PUBLIC ${MYSQL_LIBRARIES} contrib_CPP_header util Threads::Threads)
        add_test( ${test_name} ${test_name} )
    endforeach()
endif()