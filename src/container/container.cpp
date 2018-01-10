#include "container/container.h"
using namespace libMABS;

void exportContainer()
{
    // container is an abstract class and should never be initialized
    boost::python::class_<Container, std::shared_ptr<Container>>(
            "Container",
            boost::python::no_init
        )
        .def(
                "get_type_info", 
                &Container::getType
            )
        .def(
                "get_name", 
                &Container::getTypeName
            );

    //export the Nil class
    boost::python::class_<
            Nil, 
            boost::python::bases<Container>,
            std::shared_ptr<Nil>
        >("Nil")
        ;
    //tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible<
        std::shared_ptr<Nil>,
        std::shared_ptr<Container>
    >();
    // @todo @fixme python initialized the vector with one element instead of an empty vector
    boost::python::class_<ContainerVector, 
            boost::python::bases<Container>,
            std::shared_ptr<ContainerVector>
        >("ContainerVector", boost::python::init<std::shared_ptr<Container>>())
        .def(boost::python::vector_indexing_suite<
                ContainerVector,
                /*
                *    true = noproxy this means that the content of the vector is already exposed by
                *    boost python. 
                *    if this is kept as false, Container would be exposed a second time.
                *    the two Containers would be different and not inter castable.
                */
                true
            >());

    //DEPRECATED
    //make vectors of container-pointers a thing
    //IterableConverter()
    //    .from_python<std::shared_ptr<ContainerVector>>();

    //tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible< 
        std::shared_ptr<ContainerVector>,
        std::shared_ptr<Container> 
    >();
}//function