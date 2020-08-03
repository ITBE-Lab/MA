cmake_minimum_required (VERSION 3.8)

project("MS Lib" VERSION 1.1.1 DESCRIPTION "Modular Structure library")


if( PythonLibs_FOUND )
    message( "--libMS: Building shared library with python support via local pybind11." )
else()
    message( "--libMS: Building shared library without python support." )
endif()

# and for each of your dependent executable/library targets:
if( ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0) )
    # std::filesystem is still experimental and for linkage we need separated library
    message( "--libMS: Filesystem support for GCC below version 8" )
endif()

function( ms_add_module_helper target_name )
    # Get all .cpp files in the src folder 
    file(GLOB_RECURSE SRC
        "src/*.cpp"
    ) # file

    # Define library MA as target
    if( PythonLibs_FOUND )
        pybind11_add_module( "lib${target_name}" SHARED SYSTEM ${SRC} )
        target_compile_definitions( "lib${target_name}" PRIVATE WITH_PYTHON )
    else()
        # Build the plain shared library without Python support
        add_library( "lib${target_name}" SHARED ${SRC} )
    endif()

    # Activate __declspec(dllexport) for DLL export 
    # PRIVATE => depending projects will not inherit the EXPORT definition
    message( "--lib${target_name}: setting DLL_PORT_${target_name} to dllexport" )
    target_compile_definitions( "lib${target_name}" PRIVATE DLL_PORT_${target_name}=dllexport )

    # Be verbose with sql queries (in debug mode)
    #target_compile_definitions( MA PUBLIC SQL_VERBOSE )

    # CMake will add the -fPIC for compilation with GCC
    set_property(TARGET "lib${target_name}" PROPERTY POSITION_INDEPENDENT_CODE ON)

    target_link_libraries ( "lib${target_name}" PUBLIC util contrib_CPP_header )
    target_compile_definitions( "lib${target_name}" PUBLIC DLL_PORT_util=dllimport )

    # Set includes and library dependencies 
    target_include_directories("lib${target_name}" PUBLIC inc)
    
    # build the shared library in the folder of the .exe 
    set_target_properties( "lib${target_name}" PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )
    set_target_properties( "lib${target_name}" PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

    # and for each of your dependent executable/library targets:
    if( ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0) )
        # std::filesystem is still experimental and for linkage we need separated library
        target_link_libraries("lib${target_name}" PUBLIC stdc++fs )
    endif()

    # additionally link to python if available
    if( PythonLibs_FOUND )
        target_link_libraries ( "lib${target_name}" PUBLIC ${PYTHON_LIBRARIES} )
        # if we are under windows and we want to support python we need to do two things:
        if( MSVC )
            # change the name of the library (under linux lib is used as prefix automatically)
            set_target_properties( "lib${target_name}" PROPERTIES OUTPUT_NAME "lib${target_name}" )
            # change the suffix to .pdy so that python recognized the file as module 
            # (this does not effect the executable; it can still load "lib"${target_name}.pyd as dynamic linked library)
            set_target_properties( "lib${target_name}" PROPERTIES SUFFIX ".pyd" )
        else( )
            # change the name of the library (under linux lib is used as prefix automatically)
            set_target_properties( "lib${target_name}" PROPERTIES OUTPUT_NAME "${target_name}" )
        endif()

        # create the python module
        make_directory("${CMAKE_BINARY_DIR}/${target_name}")
        file(GLOB py_folder "python" )
        make_directory("${CMAKE_BINARY_DIR}/${target_name}")
        add_custom_command(
                OUTPUT "${CMAKE_BINARY_DIR}/${target_name}/copy_directory"
                COMMAND ${CMAKE_COMMAND} -E copy_directory ${py_folder} "${CMAKE_BINARY_DIR}/${target_name}"
                DEPENDS ${py_folder}
        )
        configure_file(
            ${_lib_init_py_in}
            "${CMAKE_BINARY_DIR}/generated/${target_name}/_lib_init_.py"
        )
        add_custom_command(
                OUTPUT "${CMAKE_BINARY_DIR}/${target_name}/_lib_init.py"
                COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_BINARY_DIR}/generated/${target_name}/_lib_init_.py" > 
                    "${CMAKE_BINARY_DIR}/${target_name}/_lib_init.py"
                DEPENDS "lib${target_name}" "${CMAKE_BINARY_DIR}/generated/${target_name}/_lib_init_.py"
            )
        add_custom_target(
                ${target_name} ALL
                DEPENDS "${CMAKE_BINARY_DIR}/${target_name}/_lib_init.py" 
                        "${CMAKE_BINARY_DIR}/${target_name}/copy_directory"
            )
    endif()

    # Integration of test-code, if requested
    if(BUILD_TESTING)
        # add tests for all cpp files in the tests folder
        # Get all .cpp files in the test folder
        file(GLOB_RECURSE CPP_TEST_SRC
            "tests/*.cpp"
        ) # file
        foreach(file ${CPP_TEST_SRC})
            get_filename_component(test_name ${file} NAME)
            message( "--lib${target_name}: generating test " "${target_name}-${test_name}" )
            add_executable( "${target_name}-${test_name}" ${file} )
            target_link_libraries( "${target_name}-${test_name}" PUBLIC "lib${target_name}" Threads::Threads )
            target_compile_definitions( "${target_name}-${test_name}" PUBLIC "DLL_PORT_${target_name}=dllimport" )
            add_test( "${target_name}-${test_name}" "${target_name}-${test_name}" )
        endforeach()

        # add tests for all python files in the tests folder
        if( PythonLibs_FOUND )
            find_package(PythonInterp ${requested_python_version})
            if(PythonInterp_FOUND)
                # Get all .py files in the test folder
                file(GLOB_RECURSE PY_TEST_SRC
                    "tests/*.py"
                ) # file
                foreach(file ${PY_TEST_SRC})
                    get_filename_component(test_name ${file} NAME)
                    message( "--lib${target_name}: generating test " "${target_name}-${test_name}" )
                    add_test(NAME "${target_name}-${test_name}" COMMAND ${PYTHON_EXECUTABLE} ${file} )
                endforeach()
            else()
                message("WARNING failed to generate python tests because no interpreter was found")
            endif()
        endif()
    endif()
endfunction()

file(GLOB_RECURSE _lib_init_py_in "_lib_init.py.in" )
set( _lib_init_py_in ${_lib_init_py_in} CACHE INTERNAL "" )
ms_add_module_helper( MS )
target_include_directories( libMS PUBLIC dbSQL )
target_link_libraries ( libMS PUBLIC dbSQL )

# Generate a version of MS based on the current commit

# Get the latest abbreviated commit hash of the working branch
execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    RESULT_VARIABLE RESULT_COMMIT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)

# Check if there are any modified files
execute_process(
    COMMAND git status -s
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_STATUS
    RESULT_VARIABLE RESULT_STATUS
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)

if( "${GIT_STATUS}" STREQUAL "" )
    set(GIT_STATUS "")
else()
    set(GIT_STATUS "-D")
endif()

if(NOT RESULT_COMMIT STREQUAL "0" OR NOT RESULT_STATUS STREQUAL "0")
    message( "--libMS: WARNING cannot infer verison from git (maybe git.exe is not in your environment path?)" )
    set(MS_VERSION "Unknown")
else()
    set(MS_VERSION "${PROJECT_VERSION}-${GIT_COMMIT_HASH}${GIT_STATUS}")
    message( "--libMS: Version = ${MS_VERSION}" )
endif()

configure_file( # @todo so far this only triggers during configure not at build time
    inc/ms/util/version.h.in
    ${CMAKE_BINARY_DIR}/generated/ms/util/version.h
)
target_include_directories( libMS PUBLIC ${CMAKE_BINARY_DIR}/generated)

# the function other modules shall use to create a compatitable library
function( ms_add_module target_name )
    ms_add_module_helper( ${target_name} )
    target_link_libraries ( "lib${target_name}" LINK_PUBLIC libMS )
    target_compile_definitions( "lib${target_name}" PUBLIC DLL_PORT_MS=dllimport )
endfunction()

