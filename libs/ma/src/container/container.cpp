/**
 * @file container.cpp
 * @author Markus Schmidt
 */


#include "container/container.h"
#include "util/pybind11.h"
using namespace libMA;


#ifdef WITH_PYTHON
#ifdef BOOST_PYTHON
void exportContainer( )
{
    // container is an abstract class and should never be initialized
    boost::python::class_<Container, std::shared_ptr<Container>>( "Container", boost::python::no_init );

    boost::python::class_<PyContainerVector, boost::noncopyable, boost::python::bases<Container>,
                          std::shared_ptr<PyContainerVector>>( "ContainerVector" )
        /*
         * true = noproxy this means that the content of
         * the vector is already exposed by boost python.
         * if this is kept as false, Container would be
         * exposed a second time. the two Containers would
         * be different and not inter castable.
         */
        .def( boost::python::vector_indexing_suite<PyContainerVector, true>( ) );

    // tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible<std::shared_ptr<PyContainerVector>, std::shared_ptr<Container>>( );

} // function
#else

void exportContainer( py::module& rxPyModuleId )
{
    py::class_<Container, std::shared_ptr<Container>>( rxPyModuleId, "Container" );

    // container is an abstract class and should never be initialized

    py::bind_vector_ext<PyContainerVector, Container, std::shared_ptr<PyContainerVector>>(
        rxPyModuleId, "ContainerVector", "docstr" );

    // tell boost python that pointers of these classes can be converted implicitly 
    // @todo is this still necessary?
    py::implicitly_convertible<PyContainerVector, Container>( );

} // function
#endif
#endif