cmake_minimum_required (VERSION 3.8)

project("MA Lib" VERSION 1.1.1 DESCRIPTION "The Modular Aligner Library")

ms_add_module( MA )

# Set includes and library dependencies 
target_link_libraries ( MA LINK_PUBLIC kswcpp )
if( ZLIB_FOUND )
    target_link_libraries ( MA LINK_PUBLIC minimizerscpp )
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
