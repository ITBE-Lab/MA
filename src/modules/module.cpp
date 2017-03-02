#include "modules/module.h"

std::shared_ptr<Container> Module::execute(std::shared_ptr<Container> pInput) 
{
    std::cout << "i should not be called" << std::endl; 
    return pInput;
}


std::shared_ptr<Container> Printer::execute(std::shared_ptr<Container> pInput) 
{
    std::cout << "printer" << std::endl; 
    return pInput;
}

void exportModule()
{
	boost::python::class_<Module>("Module");
	boost::python::class_<Printer, boost::python::bases<Module>>("Printer");
}