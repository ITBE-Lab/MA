cmake_minimum_required (VERSION 3.8)

project("MA Lib" VERSION 1.1.1 DESCRIPTION "The Modular Aligner Library")

ms_add_module( MA )

# Set includes and library dependencies 
target_link_libraries ( MA LINK_PUBLIC kswcpp contrib_CPP_header )
if( ZLIB_FOUND )
    target_link_libraries ( MA LINK_PUBLIC minimizerscpp )
endif()

# and for each of your dependent executable/library targets:
if( ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0) )
    # std::filesystem is still experimental and for linkage we need separated library
    message( "libMA: Filesystem support for GCC below version 8" )
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
    message( "--libMA: Building shared library components that rely on zlib." )
else()
    message( "--libMA: Not building shared library components that rely on zlib." )
endif()
