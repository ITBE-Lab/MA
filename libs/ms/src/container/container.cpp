/**
 * @file container.cpp
 * @author Markus Schmidt
 */


#include "ms/container/container.h"
#include "ms/util/pybind11.h"

using namespace libMS;


#ifdef WITH_PYTHON

void exportContainer( SubmoduleOrganizer& xOrganizer )
{
    py::class_<Container, std::shared_ptr<Container>>( xOrganizer.container( ), "Container" ).def( py::init<>( ) );

    // container is an abstract class and should never be initialized

    py::bind_vector_ext<PyContainerVector, Container, std::shared_ptr<PyContainerVector>>(
        xOrganizer.container( ), "ContainerVector", "docstr" );

    // tell boost python that pointers of these classes can be converted implicitly
    py::implicitly_convertible<PyContainerVector, Container>( );

} // function
#endif