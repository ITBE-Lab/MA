/**
 * @file export.h
 * @brief defines a class to organize the python module.
 */
#pragma once

#ifdef WITH_PYTHON
// Bug in Python 3.7 and 3.8 that breaks nlohmann::json
// See:(https://bugs.python.org/issue36020)
#if defined( WIN32 ) && !defined( HAVE_SNPRINTF ) && defined( _MSC_VER ) && _MSC_VER >= 1900
#define HAVE_SNPRINTF
#endif

#include <pybind11/pybind11.h>
namespace py = pybind11;

namespace libMS
{

#if defined( __GNUC__ )
#pragma GCC visibility push( hidden ) // @todo is there a way to resolve this warning?
#endif

/**
 * @brief creates several submodules in python.
 * @details
 * Some of these submodules are then imported by _lib_init.py.
 * Submodule organisation:
 * - modules: imported; contains all modules; Is used automatically by exportModule; Do not use this for anything bu
 *            modules.
 * - _modules: not imported; contains the raw cpp modules (not wrapped by the python to cpp interface).
 * - containers: imported; use this to expose any containers that shall be visible to python.
 * - _containers: not imported; use this to expose any containers that shall not be visible to python (g.e. return 
 *                types of functions that are not initiable from python directly).
 * - util: imported; any type of utility functions / classes
 * - _util: not imported; any type of hidden utility functions / classes (e.g. hidden return types of functions).
 */
class SubmoduleOrganizer
{
  private:
    py::module xHiddenModules;
    py::module xHiddenContainer;
    py::module xHiddenUtil;
    // order matters
    py::module xModules;
    py::module xContianers;
    py::module xUtility;

  public:
    SubmoduleOrganizer( py::module& rMainModule )
        : xHiddenModules( rMainModule.def_submodule( "_modules", "Private Modules" ) ),
          xHiddenContainer( rMainModule.def_submodule( "_containers", "Private Containers" ) ),
          xHiddenUtil( rMainModule.def_submodule( "_util", "Private Utility" ) ),
          xModules( rMainModule.def_submodule( "modules", "Modules that implement various algorithms" ) ),
          xContianers( rMainModule.def_submodule( "containers", "Containers that hold various datatypes" ) ),
          xUtility( rMainModule.def_submodule( "util", "Utility functions and classes" ) )
    {} // constructor

    inline py::module& module( )
    {
        return xModules;
    } // method


    inline py::module& container( )
    {
        return xContianers;
    } // method

    inline py::module& util( )
    {
        return xUtility;
    } // method

    inline py::module& _module( )
    {
        return xHiddenModules;
    } // method

    inline py::module& _container( )
    {
        return xHiddenContainer;
    } // method

    inline py::module& _util( )
    {
        return xHiddenUtil;
    } // method

}; // class

#if defined( __GNUC__ )
#pragma GCC visibility pop
#endif

} // namespace libMS

#endif