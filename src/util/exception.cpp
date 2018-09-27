/**
 * @file exception.cpp
 * @author Markus Schmidt
 */
#include "util/exception.h"

#ifdef WITH_PYTHON
void translator( ModuleIO_Exception const &x )
{
    PyErr_SetString( PyExc_RuntimeError, x.what( ) );
}

void exportExceptions( )
{
    boost::python::register_exception_translator<ModuleIO_Exception>( translator );
} // function
#endif