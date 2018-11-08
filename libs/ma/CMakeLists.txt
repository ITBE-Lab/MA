cmake_minimum_required (VERSION 3.8)

project("MA Lib" VERSION 1.0.1 DESCRIPTION "Core library of MA")

# Look for Python
find_package( PythonLibs 3.5 )
if( PythonLibs_FOUND )
	message( "Found Python Libs")
endif()

# Get all .cpp files in the src folder
file(GLOB MA_SRC "src/*/*.cpp" )

# Define library MA as target
if( PythonLibs_FOUND )
    pybind11_add_module( MA SHARED SYSTEM ${MA_SRC} )
    target_compile_definitions( MA PRIVATE WITH_PYTHON )
else()
    add_library( MA SHARED ${MA_SRC} ) 
endif()

# Activate __declspec(dllexport) for DLL export 
# PRIVATE => depending projects will not inherit the EXPORT definition
target_compile_definitions( MA PRIVATE EXPORT )

# CMake will add the -fPIC for compilation with GCC
set_property(TARGET MA PROPERTY POSITION_INDEPENDENT_CODE ON)

# Set includes and library dependencies 
target_include_directories( MA PUBLIC inc )
target_link_libraries ( MA PUBLIC kswcpp CppSQLite3 )

# additionally link to python if available
if( PythonLibs_FOUND )
    target_link_libraries ( MA PUBLIC ${PYTHON_LIBRARIES} )
endif()