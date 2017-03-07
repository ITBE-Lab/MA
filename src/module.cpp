#include "module.h"

void exportModule()
{
    //module is an abstract class and should never be initialized
	boost::python::class_<Module, std::shared_ptr<Module>>("Module", boost::python::no_init)
        .def("execute", &Module::saveExecute)
        .def("getInputType", &Module::getInputType)
        .def("getOutputType", &Module::getOutputType)
    ;
    //test class
	boost::python::class_<Printer, boost::python::bases<Module>, std::shared_ptr<Printer>>("Printer");

    //tell boost python that it's possible to convert shared pointers with these classes
    boost::python::implicitly_convertible<std::shared_ptr<Printer>,std::shared_ptr<Module>>();

}