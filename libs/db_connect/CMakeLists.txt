cmake_minimum_required (VERSION 3.8)

project("db_connect" VERSION 0.1 DESCRIPTION "Database connection library")

INCLUDE(cmake/FindMySQL.cmake)

# Check for PostgreSQL Support
# (With Windows, PostgreSQL DLL location must be added to search path. (bin and lib folders of installation)) 
find_package(PostgreSQL)

# Define the target dbSQL
add_library(dbSQL INTERFACE)
target_compile_definitions( dbSQL INTERFACE DLL_PORT_util=dllimport )

# link to the fort object file
target_link_libraries ( dbSQL INTERFACE Threads::Threads fort util contrib_CPP_header )
# location of all includes 
target_include_directories( dbSQL INTERFACE connectors . )

if( MySQL_FOUND )
    message( "Found MySQL includes and libraries")
    target_link_libraries ( dbSQL INTERFACE ${MYSQL_LIBRARIES} )
    target_include_directories( dbSQL INTERFACE ${MYSQL_INCLUDE_DIRS} )
endif()

if( PostgreSQL_FOUND )
    message( "Found PostgreSQL - library: " ${PostgreSQL_LIBRARIES} )
    target_link_libraries ( dbSQL INTERFACE ${PostgreSQL_LIBRARIES} )
    target_include_directories( dbSQL INTERFACE ${PostgreSQL_INCLUDE_DIRS} )
endif()

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
        target_link_libraries( ${test_name} PUBLIC dbSQL)
		# and for each of your dependent executable/library targets:
		if( ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0) )
		    # std::filesystem is still experimental and for linkage we need separated library
		    message( "-- libMA: Filesystem support for GCC below version 8" )
		    target_link_libraries( ${test_name} PUBLIC stdc++fs )
		endif()
        add_test( ${test_name} ${test_name} )
    endforeach()

    if( PostgreSQL_FOUND )
        # Add all PostgreSQL test
        file(GLOB POSTGRESQL_TEST_SRC
            "tests/PostgreSQL/*.cpp"
        ) # file
        
        foreach(file ${POSTGRESQL_TEST_SRC})
            get_filename_component(test_name ${file} NAME)
            message( "-- Add PostgreSQL test: " ${test_name} )
            add_executable( ${test_name} ${file} )
	    	# Add include directories to your target. PRIVATE is useful with multi-target projects - see documentation of target_include_directories for more info
            target_include_directories( ${test_name} PUBLIC dbSQL connectors )

            # Add libraries to link your target against. Again, PRIVATE is important for multi-target projects
            target_link_libraries( ${test_name} PUBLIC dbSQL )
            # and for each of your dependent executable/library targets:
		    if( ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0) )
		        # std::filesystem is still experimental and for linkage we need separated library
		        message( "-- libMA: Filesystem support for GCC below version 8" )
		        target_link_libraries( ${test_name} PUBLIC stdc++fs )
		    endif()
            add_test( ${test_name} ${test_name} )
        endforeach()
    endif()
endif()


