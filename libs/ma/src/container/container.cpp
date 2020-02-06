/**
 * @file container.cpp
 * @author Markus Schmidt
 */


#include "container/container.h"
#include "util/pybind11.h"
using namespace libMA;


#ifdef WITH_PYTHON

void exportContainer( py::module& rxPyModuleId )
{
    py::class_<Container, std::shared_ptr<Container>>( rxPyModuleId, "Container" );

    // container is an abstract class and should never be initialized

    py::bind_vector_ext<PyContainerVector, Container, std::shared_ptr<PyContainerVector>>(
        rxPyModuleId, "ContainerVector", "docstr" );

    // tell boost python that pointers of these classes can be converted implicitly 
    py::implicitly_convertible<PyContainerVector, Container>( );

} // function
#endif