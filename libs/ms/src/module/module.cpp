/**
 * @file module.cpp
 * @author Markus Schmidt
 */
#include "module/module.h"
using namespace libMS;

#ifdef WITH_PYTHON

void exportModuleClass( py::module& rxPyModuleId )
{
    /*
     * Module and VolatileModule can be extended by python and execute can be overloaded
     * HOWEVER: the default constructor of Module and VolatileModule MUST be calles otherwise the programm segfaults
     */
    py::class_<PyModule<false>, ModuleWrapperPyToCpp<false>, std::shared_ptr<PyModule<false>>>( rxPyModuleId, "Module" )
        .def( py::init<>( ) ) // default constructor
        .def( "execute", &PyModule<false>::execute );

    py::class_<PyModule<true>, ModuleWrapperPyToCpp<true>, std::shared_ptr<PyModule<true>>>( rxPyModuleId,
                                                                                             "VolatileModule" )
        .def( py::init<>( ) ) // default constructor
        .def( "execute", &PyModule<true>::execute );

    py::class_<BasePledge, std::shared_ptr<BasePledge>>( rxPyModuleId, "BasePledge" );

    py::class_<PyPledgeVector, std::shared_ptr<PyPledgeVector>>( rxPyModuleId, "VectorPledge" )
        .def( py::init<>( ) ) // default constructor
        .def( "append", &PyPledgeVector::append )
        .def( "get", &PyPledgeVector::get )
        .def( "clear", &PyPledgeVector::clear )
        .def( "simultaneous_get", &PyPledgeVector::simultaneousGetPy );

    py::implicitly_convertible<PyPledgeVector, BasePledge>( );


    typedef Pledge<Container, false, PyPledgeVector> TP_MODULE_PLEDGE;
    py::class_<TP_MODULE_PLEDGE, BasePledge, std::shared_ptr<TP_MODULE_PLEDGE>>( rxPyModuleId, "ModulePledge" )
        .def( py::init<std::shared_ptr<PyModule<false>>, std::shared_ptr<PyPledgeVector>>( ) )
        .def("exec_time", &TP_MODULE_PLEDGE::execTime)
        .def( "get", &TP_MODULE_PLEDGE::get );

    py::implicitly_convertible<TP_MODULE_PLEDGE, BasePledge>( );


    typedef Pledge<Container, true, PyPledgeVector> TP_VOLATILE_PLEDGE;
    py::class_<TP_VOLATILE_PLEDGE, BasePledge, std::shared_ptr<TP_VOLATILE_PLEDGE>>( rxPyModuleId,
                                                                                     "VolatileModulePledge" )
        .def( py::init<std::shared_ptr<PyModule<true>>, std::shared_ptr<PyPledgeVector>>( ) )
        .def("exec_time", &TP_VOLATILE_PLEDGE::execTime)
        .def( "get", &TP_VOLATILE_PLEDGE::get );

    py::implicitly_convertible<TP_VOLATILE_PLEDGE, BasePledge>( );


    typedef Pledge<Container, false> TP_PLEDGE;
    py::class_<TP_PLEDGE, BasePledge, std::shared_ptr<TP_PLEDGE>>( rxPyModuleId, "Pledge" )
        .def( py::init<>( ) ) // default constructor
        .def("exec_time", &TP_PLEDGE::execTime)
        .def("wait_on_lock_time", &TP_PLEDGE::waitTime)
        .def( "set", &TP_PLEDGE::set )
        .def( "get", &TP_PLEDGE::get );
    py::implicitly_convertible<TP_PLEDGE, BasePledge>( );
} // function

#endif