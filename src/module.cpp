#include "module.h"

void exportModule()
{
    //module is an abstract class and should never be initialized
	boost::python::class_<Module>("Module", boost::python::no_init)
        .def("execute", &Module::saveExecute)
        .def("getInputType", &Module::getInputType)
        .def("getOutputType", &Module::getOutputType)
    ;
    //test class
	boost::python::class_<Printer, boost::python::bases<Module>>("Printer");

}