/**
 * @file exception.cpp
 * @author Markus Schmidt
 */
#include "util/exception.h"

#ifdef WITH_PYTHON
#ifdef BOOST_PYTHON
void translator( AnnotatedException const& x )
{
    PyErr_SetString( PyExc_RuntimeError, x.what( ) );
}

void exportExceptions( )
{
    boost::python::register_exception_translator<AnnotatedException>( translator );
} // function
#endif
#endif