/**
 * @file exception.cpp
 * @author Markus Schmidt
 */
#include "util/exception.h"

#ifdef WITH_PYTHON
void translator( AnnotatedException const& x )
{
    PyErr_SetString( PyExc_RuntimeError, x.what( ) );
}
#ifdef BOOST_PYTHON

void exportExceptions( )
{
    boost::python::register_exception_translator<AnnotatedException>( translator );
} // function
#else
void exportExceptions( py::module& rxPyModuleId )
{
    py::register_exception<AnnotatedException>( rxPyModuleId, "AnnotatedException" );
} // function
#endif
#endif