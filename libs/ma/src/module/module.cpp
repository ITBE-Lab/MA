/**
 * @file module.cpp
 * @author Markus Schmidt
 */
#include "module/module.h"
#include "module/splitter.h"
using namespace libMA;

#ifdef WITH_PYTHON
#ifdef BOOST_PYTHON
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
        //.def( "set", &PyPledgeVector::set )
        .def( "append", &PyPledgeVector::append )
        .def( "get", &PyPledgeVector::get )
        .def( "simultaneous_get", &PyPledgeVector::simultaneousGetPy )
        //.def_readwrite( "exec_time", &PyPledgeVector::execTime )
        ;
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
#else

void exportModuleClass( py::module& rxPyModuleId )
{
    // module is an abstract class and should never be initialized
    py::class_<PyModule<false>, std::shared_ptr<PyModule<false>>>( rxPyModuleId, "Module" )
        .def( "execute", &PyModule<false>::execute );
    py::class_<PyModule<true>, std::shared_ptr<PyModule<true>>>( rxPyModuleId, "VolatileModule" )
        .def( "execute", &PyModule<true>::execute );

    py::class_<BasePledge, std::shared_ptr<BasePledge>>( rxPyModuleId, "BasePledge" );

    py::class_<PyPledgeVector, std::shared_ptr<PyPledgeVector>>( rxPyModuleId, "VectorPledge" )
        .def( py::init<>( ) ) // default constructor
        .def( "append", &PyPledgeVector::append )
        .def( "get", &PyPledgeVector::get )
        .def( "simultaneous_get", &PyPledgeVector::simultaneousGetPy );

    py::implicitly_convertible<PyPledgeVector, BasePledge>( );


    typedef Pledge<Container, false, PyPledgeVector> TP_MODULE_PLEDGE;
    py::class_<TP_MODULE_PLEDGE, BasePledge, std::shared_ptr<TP_MODULE_PLEDGE>>( rxPyModuleId, "ModulePledge" )
        .def( py::init<std::shared_ptr<PyModule<false>>, std::shared_ptr<PyPledgeVector>>( ) )
        .def( "get", &TP_MODULE_PLEDGE::get );

    py::implicitly_convertible<TP_MODULE_PLEDGE, BasePledge>( );


    typedef Pledge<Container, true, PyPledgeVector> TP_VOLATILE_PLEDGE;
    py::class_<TP_VOLATILE_PLEDGE, BasePledge, std::shared_ptr<TP_VOLATILE_PLEDGE>>( rxPyModuleId, 
                                                                                     "VolatileModulePledge" )
        .def( py::init<std::shared_ptr<PyModule<true>>, std::shared_ptr<PyPledgeVector>>( ) );

    py::implicitly_convertible<TP_VOLATILE_PLEDGE, BasePledge>( );


    typedef Pledge<Container, false> TP_PLEDGE;
    py::class_<TP_PLEDGE, std::shared_ptr<TP_PLEDGE>>( rxPyModuleId, "Pledge" )
        .def( py::init<>( ) ) // default constructor
        .def( "set", &TP_PLEDGE::set )
        .def( "get", &TP_PLEDGE::get );
    py::implicitly_convertible<TP_PLEDGE, BasePledge>( );
} // function

#endif
#endif