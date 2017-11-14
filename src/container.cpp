#include "container.h"



void exportContainer()
{
    // container is an abstract class and should never be initialized
	boost::python::class_<Container, std::shared_ptr<Container>>(
            "Container", 
            "abstract\n"
            "Any class holding data should inherit Container.\n",
            boost::python::no_init
        )
        .def(
                "get_type_info", 
                &Container::getType,
                "arg1: self\n"
                "returns: an enum describing the type of data stored.\n"
                "\n"
                "Used to check weather a module can work with the given containers.\n"
            );

            
    boost::python::class_<std::vector<std::shared_ptr<Container>>>(
            "ContainerVector"
        )
        .def(boost::python::vector_indexing_suite<
                std::vector<std::shared_ptr<Container>>,
                /*
                *	true = noproxy this means that the content of the vector is already exposed by
                *	boost python. 
                *	if this is kept as false, SeedsVector would be exposed a second time.
                *	the two SeedsVector would be different and not inter castable.
                */
                true
            >());

    //make vectors of container-pointers a thing
    iterable_converter()
        .from_python<std::vector<std::shared_ptr<Container>>>();

}//function