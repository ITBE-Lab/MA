/**
 * @file module.cpp
 * @author Markus Schmidt
 */
#include "module/module.h"
#include "module/splitter.h"
using namespace libMA;

#ifdef WITH_PYTHON
void exportModule( )
{
    // module is an abstract class and should never be initialized
    boost::python::class_<Module>( "Module", boost::python::no_init )
        .def( "execute", &Module::pyExecute )
        .def( "get_input_type", &Module::getInputType )
        .def( "get_output_type", &Module::getOutputType )
        .def( "get_name", &Module::getName )
        .def( "promise_me", &Module::promiseMe
              //,boost::python::with_custodian_and_ward_postcall<0,1,
              //    boost::python::with_custodian_and_ward_postcall<0,2>
              //>()
        );
    boost::python::class_<Pledge, boost::noncopyable, boost::python::bases<Container>,
                          std::shared_ptr<Pledge>>(
        "Pledge", boost::python::init<std::shared_ptr<Container>>( ) )
        .def( "make_pledge", &Pledge::makePyPledge
              //,boost::python::with_custodian_and_ward_postcall<0,1,
              //    boost::python::with_custodian_and_ward_postcall<0,3>
              //>()
              )
        .staticmethod( "make_pledge" )
        .def( "set", &Pledge::set_boost )
        .def( "get", &Pledge::get
              //,boost::python::with_custodian_and_ward_postcall<1,0>()
              )
        .def( "clear_graph", &Pledge::clear_graph
              //,boost::python::with_custodian_and_ward_postcall<1,0>()
              )
        .def( "simultaneous_get", &Pledge::simultaneousGet )
        .def( "simultaneous_get", &simGetBoost1 )
        .def( "simultaneous_get", &simGetBoost2 )
        .staticmethod( "simultaneous_get" )
        .def( "get_pledger", &Pledge::getPledger )
        .def( "get_graph_desc", &Pledge::getGraphDesc )
        .def_readwrite( "exec_time", &Pledge::execTime );

    // tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible<std::shared_ptr<Pledge>, std::shared_ptr<Container>>( );

    IterableConverter( ).from_python<std::vector<std::shared_ptr<Pledge>>>( );
} // function
#endif