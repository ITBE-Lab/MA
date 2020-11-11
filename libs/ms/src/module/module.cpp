/**
 * @file module.cpp
 * @author Markus Schmidt
 */
#define USE_DLL_EXPORT
#include "ms/module/module.h"
using namespace libMS;

#ifdef _MSC_VER
__declspec(dllexport) const size_t BasePledge::uiDefaultGraphThread = 0;
__declspec(dllexport) size_t BasePledge::uiThreadCurrentlyBuildingGraph = uiDefaultGraphThread;
#else
const size_t BasePledge::uiDefaultGraphThread = 0;
size_t BasePledge::uiThreadCurrentlyBuildingGraph = uiDefaultGraphThread;
#endif

#ifdef WITH_PYTHON

void exportModuleClass( SubmoduleOrganizer& xOrganizer )
{
    /*
     * Module and VolatileModule can be extended by python and execute can be overloaded
     * HOWEVER: the default constructor of Module and VolatileModule MUST be calles otherwise the programm segfaults
     */
    py::class_<PyModule<false>, ModuleWrapperPyToCpp<false>, std::shared_ptr<PyModule<false>>>( xOrganizer.util( ),
                                                                                                "Module" )
        .def( py::init<>( ) ) // default constructor
        .def( "execute", &PyModule<false>::execute );

    py::class_<PyModule<true>, ModuleWrapperPyToCpp<true>, std::shared_ptr<PyModule<true>>>( xOrganizer.util( ),
                                                                                             "VolatileModule" )
        .def( py::init<>( ) ) // default constructor
        .def( "execute", &PyModule<true>::execute );

    py::class_<BasePledge, std::shared_ptr<BasePledge>>( xOrganizer.util( ), "BasePledge" )
        .def_readonly_static( "default_graph_thread", &BasePledge::uiDefaultGraphThread )
        .def_readwrite_static( "current_graph_thread", &BasePledge::uiThreadCurrentlyBuildingGraph );

    py::class_<PyPledgeVector, std::shared_ptr<PyPledgeVector>>( xOrganizer.util( ), "VectorPledge" )
        .def( py::init<>( ) ) // default constructor
        .def( "append", &PyPledgeVector::append )
        .def( "get", &PyPledgeVector::get )
        .def( "clear", &PyPledgeVector::clear )
        .def( "simultaneous_get", &PyPledgeVector::simultaneousGetPy );

    py::implicitly_convertible<PyPledgeVector, BasePledge>( );


    typedef Pledge<Container, false, PyPledgeVector> TP_MODULE_PLEDGE;
    py::class_<TP_MODULE_PLEDGE, BasePledge, std::shared_ptr<TP_MODULE_PLEDGE>>( xOrganizer._util( ), "ModulePledge" )
        .def( py::init<std::shared_ptr<PyModule<false>>, std::shared_ptr<PyPledgeVector>>( ) )
        .def( "exec_time", &TP_MODULE_PLEDGE::execTime )
        .def( "get", &TP_MODULE_PLEDGE::get );

    py::implicitly_convertible<TP_MODULE_PLEDGE, BasePledge>( );


    typedef Pledge<Container, true, PyPledgeVector> TP_VOLATILE_PLEDGE;
    py::class_<TP_VOLATILE_PLEDGE, BasePledge, std::shared_ptr<TP_VOLATILE_PLEDGE>>( xOrganizer._util( ),
                                                                                     "VolatileModulePledge" )
        .def( py::init<std::shared_ptr<PyModule<true>>, std::shared_ptr<PyPledgeVector>>( ) )
        .def( "exec_time", &TP_VOLATILE_PLEDGE::execTime )
        .def( "get", &TP_VOLATILE_PLEDGE::get );

    py::implicitly_convertible<TP_VOLATILE_PLEDGE, BasePledge>( );


    typedef Pledge<Container, false> TP_PLEDGE;
    py::class_<TP_PLEDGE, BasePledge, std::shared_ptr<TP_PLEDGE>>( xOrganizer.util( ), "Pledge" )
        .def( py::init<>( ) ) // default constructor
        .def( "exec_time", &TP_PLEDGE::execTime )
        .def( "wait_on_lock_time", &TP_PLEDGE::waitTime )
        .def( "set", &TP_PLEDGE::set )
        .def( "get", &TP_PLEDGE::get );
    py::implicitly_convertible<TP_PLEDGE, BasePledge>( );

#if 0 // @todo needs it's own class that holds a python object
    typedef Pledge<Container, false, PyContainerVector> TP_PYTHON_PLEDGE;
    py::class_<TP_PYTHON_PLEDGE, BasePledge, std::shared_ptr<TP_PYTHON_PLEDGE>>( xOrganizer._util( ), "PythonPledge" )
        .def( py::init<std::shared_ptr<PyModule<false>>, std::shared_ptr<PyPledgeVector>>( ) )
        .def( "exec_time", &TP_PYTHON_PLEDGE::execTime )
        .def( "get", &TP_PYTHON_PLEDGE::get );
#endif
} // function

#endif