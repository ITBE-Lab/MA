cmake_minimum_required (VERSION 3.8)

project("MA Lib" VERSION 1.0.1 DESCRIPTION "Core library of MA")

# Look for Python
find_package( PythonLibs 3.5 )
if( PythonLibs_FOUND )
	message( "Found Python Libs")
endif()

# Look for boost
# Info: https://feralchicken.wordpress.com/2013/12/07/boost-python-hello-world-example-using-cmake/
# If your boost is not within the standard path, then set BOOST_ROOT appropriately
set(BOOST_ROOT "/opt/dev/")
find_package( Boost 1.60.0 COMPONENTS python3 )

# Define library libMA as target
add_library(libMA SHARED
	src/container/alignment.cpp
	src/container/bwt_large.cpp
	src/container/container.cpp
	src/container/fMIndex.cpp
	src/container/is.cpp
	src/container/nucSeq.cpp
	src/container/pack.cpp
	src/container/qSufSort.cpp
	src/container/seed.cpp
	src/container/segment.cpp
	src/container/soc.cpp
    
	src/module/binarySeeding.cpp
	src/module/dbWriter.cpp
	src/module/fileReader.cpp
	src/module/fileWriter.cpp
	src/module/harmonization.cpp
	src/module/mappingQuality.cpp
	src/module/module.cpp
	src/module/needlemanWunsch.cpp
	src/module/otherSeeding.cpp
	src/module/pairedReads.cpp
	src/module/splitter.cpp
	src/module/stripOfConsideration.cpp
	
	src/sample_consensus/ransac.cpp
	src/sample_consensus/sac.cpp
	src/sample_consensus/sac_model.cpp
	src/sample_consensus/sac_model_line.cpp
	src/sample_consensus/sac_model_plane.cpp
	src/sample_consensus/test_ransac.cpp
	
	src/util/default_parameters.cpp
	src/util/exception.cpp
	src/util/export.cpp
	src/util/support.cpp
	src/util/system.cpp
	src/util/threadPool.cpp
) # add_library

if( MSVC )
	# Activate __declspec(dllexport) for DLL export 
	# PRIVATE => depending projects will not inherit the EXPORT definition
	target_compile_definitions( libMA PRIVATE EXPORT )  
endif()

# CMake will add the -fPIC for compilation with GCC
set_property(TARGET libMA PROPERTY POSITION_INDEPENDENT_CODE ON)

# Set includes and library dependencies 
target_include_directories(libMA PUBLIC inc ../kswcpp/inc)
target_link_libraries ( libMA LINK_PUBLIC kswcpp )

# Add the Python support to libMA (if Python and Boost are available) 
if( PythonLibs_FOUND AND Boost_FOUND )
	message( "Add Python Support to MA")
	if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
		# GCC complains that boost uses the deprecated auto_ptr.
		# Corresponding warning will be switched off.
		message( "Add -Wnodeprecated-declarations for gcc")
		target_compile_options( libMA PRIVATE "-Wno-deprecated-declarations")
	endif()
	
	target_include_directories( libMA PRIVATE ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} )
	target_compile_definitions( libMA PRIVATE BOOST_ALL_DYN_LINK WITH_PYTHON )
	target_link_libraries( libMA PUBLIC ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} )
endif()