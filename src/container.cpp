#include "container.h"

void exportContainer()
{
    //module is an abstract class and should never be initialized
	boost::python::class_<Container>("Container", boost::python::no_init);
}