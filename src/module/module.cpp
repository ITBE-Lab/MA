/**
 * @file module.cpp
 * @author Markus Schmidt
 */
#include "module/module.h"
#include "module/splitter.h"
using namespace libMA;

#ifdef WITH_PYTHON
void exportModuleClass( )
{
    // module is an abstract class and should never be initialized
    boost::python::class_<PyModule<false>>( "Module", boost::python::no_init )
        .def( "execute", &PyModule<false>::execute );
    boost::python::class_<PyModule<true>>( "VolatileModule", boost::python::no_init )
        .def( "execute", &PyModule<true>::execute );

    boost::python::class_<PyPledge, boost::noncopyable, std::shared_ptr<PyPledge>>( "Pledge" )
        .def( "set", &PyPledge::set )
        .def( "get", &PyPledge::get )
        .def( "simultaneous_get", &PyPledge::simultaneousGet )
        .staticmethod( "simultaneous_get" )
        .def_readwrite( "exec_time", &PyPledge::execTime );
} // function
#endif