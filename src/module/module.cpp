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


    boost::python::class_<BasePledge, boost::noncopyable, std::shared_ptr<BasePledge>>( "BasePledge",
                                                                                        boost::python::no_init );


    boost::python::class_<PyPledgeVector, boost::noncopyable, std::shared_ptr<PyPledgeVector>>( "VectorPledge" )
        .def( "set", &PyPledgeVector::set )
        .def( "append", &PyPledgeVector::append )
        .def( "get", &PyPledgeVector::get )
        .def( "simultaneous_get", &PyPledgeVector::simultaneousGet )
        .staticmethod( "simultaneous_get" )
        .def_readwrite( "exec_time", &PyPledgeVector::execTime );
    boost::python::implicitly_convertible<std::shared_ptr<PyPledgeVector>, std::shared_ptr<BasePledge>>( );

    typedef Pledge<Container, false, PyPledgeVector> TP_MODULE_PLEDGE;
    boost::python::class_<TP_MODULE_PLEDGE, boost::noncopyable, std::shared_ptr<TP_MODULE_PLEDGE>>(
        "ModulePledge", boost::python::init<std::shared_ptr<PyModule<false>>, std::shared_ptr<PyPledgeVector>>( ) )
            .def( "get", &TP_MODULE_PLEDGE::get );
    boost::python::implicitly_convertible<std::shared_ptr<TP_MODULE_PLEDGE>, std::shared_ptr<BasePledge>>( );

    typedef Pledge<Container, true, PyPledgeVector> TP_VOLATILE_PLEDGE;
    boost::python::class_<TP_VOLATILE_PLEDGE, boost::noncopyable, std::shared_ptr<TP_VOLATILE_PLEDGE>>(
        "VolatileModulePledge",
        boost::python::init<std::shared_ptr<PyModule<true>>, std::shared_ptr<PyPledgeVector>>( ) )
            .def( "get", &TP_VOLATILE_PLEDGE::get );
    boost::python::implicitly_convertible<std::shared_ptr<TP_VOLATILE_PLEDGE>, std::shared_ptr<BasePledge>>( );

    typedef Pledge<Container, false> TP_PLEDGE;
    boost::python::class_<TP_PLEDGE, boost::noncopyable, std::shared_ptr<TP_PLEDGE>>( "Pledge" )
        .def( "set", &TP_PLEDGE::set )
        .def( "get", &TP_PLEDGE::get );
    boost::python::implicitly_convertible<std::shared_ptr<TP_PLEDGE>, std::shared_ptr<BasePledge>>( );
} // function
#endif