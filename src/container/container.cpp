/**
 * @file container.cpp
 * @author Markus Schmidt
 */
#include "container/container.h"
using namespace libMA;


#ifdef WITH_PYTHON
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
#endif