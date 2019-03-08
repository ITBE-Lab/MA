/**
 * @file module.cpp
 * @author Markus Schmidt
 */
#include "module/module.h"
#include "module/splitter.h"
using namespace libMA;

#ifdef WITH_PYTHON

void exportModuleClass( py::module& rxPyModuleId )
{
    // module is an abstract class and should never be initialized
    py::class_<PyModule<false>, std::shared_ptr<PyModule<false>>>( rxPyModuleId, "Module" )
        .def( "execute", &PyModule<false>::execute )
        .def( "is_finished", &PyModule<false>::isFinished );
    py::class_<PyModule<true>, std::shared_ptr<PyModule<true>>>( rxPyModuleId, "VolatileModule" )
        .def( "execute", &PyModule<true>::execute )
        .def( "is_finished", &PyModule<true>::isFinished );

    py::class_<BasePledge, std::shared_ptr<BasePledge>>( rxPyModuleId, "BasePledge" )
        .def( "is_finished", &BasePledge::isFinished );

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
        .def( py::init<std::shared_ptr<PyModule<true>>, std::shared_ptr<PyPledgeVector>>( ) )
        .def( "get", &TP_VOLATILE_PLEDGE::get );

    py::implicitly_convertible<TP_VOLATILE_PLEDGE, BasePledge>( );


    typedef Pledge<Container, false> TP_PLEDGE;
    py::class_<TP_PLEDGE, BasePledge, std::shared_ptr<TP_PLEDGE>>( rxPyModuleId, "Pledge" )
        .def( py::init<>( ) ) // default constructor
        .def( "set", &TP_PLEDGE::set )
        .def( "get", &TP_PLEDGE::get );
    py::implicitly_convertible<TP_PLEDGE, BasePledge>( );
} // function

#endif