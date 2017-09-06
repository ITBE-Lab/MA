#include "module.h"

void exportModule()
{
    //module is an abstract class and should never be initialized
	boost::python::class_<Module>(
            "Module", 
            "class: Module\n"
            "   Abstract Class.\n"
            "   Any class implementing a algorithm should inherit Module\n",
            boost::python::no_init
        )
        .def(
                "execute",
                &Module::saveExecute,
                "method: execute(in)\n"
                "   in: tuple/vector of containers. "
                "Their types must match the input type of this module\n"
                "   returns: a container holding the results of the "
                "computations done by this module. Must match the output type of the module.\n"
                "\n"
                "   perform the task implemented by the respective instance.\n"
            )
        .def(
                "get_input_type",
                &Module::getInputType,
                "method: get_input_type()\n"
                "   returns: the type of containers that is expected as input by this module\n"

            )
        .def(
                "get_output_type",
                &Module::getOutputType,
                "method: get_output_type()\n"
                "   returns: the type of containers that is expected as output from this module\n"
            )
    ;
    //test class
	boost::python::class_<Printer, boost::python::bases<Module>>(
            "Printer",
            "class: Printer\n"
            "   prints the given Containers to std::out\n"
        );

}