/** 
 * @file exception.cpp
 * @author Markus Schmidt
 */
#include "util/exception.h"

void translator(ModuleIO_Exception const& x) {
    PyErr_SetString(PyExc_RuntimeError, x.what()); 
}

#ifdef WITH_PYTHON
void exportExceptions(){
     boost::python::register_exception_translator<ModuleIO_Exception>(translator);
}//function
#endif