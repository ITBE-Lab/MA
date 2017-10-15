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
        ;

    //make vectors of container-pointers a thing
    iterable_converter()
        .from_python<std::vector<std::shared_ptr<Container>>>()
        .from_python<std::vector<ContainerType>>();

    //export the containertype enum
    boost::python::enum_<ContainerType>("CppContainerType")
        .value("fM_index", ContainerType::fM_index)
        .value("nucSeq", ContainerType::nucSeq)
        .value("alignment", ContainerType::alignment)
        .value("packedNucSeq", ContainerType::packedNucSeq)
        .value("segmentList", ContainerType::segmentList)
        .value("segment", ContainerType::segment)
        .value("seed", ContainerType::seed)
        .value("seeds", ContainerType::seeds)
        .value("seedsVector", ContainerType::seedsVector)
        .value("sa_interval", ContainerType::sa_interval)
        .value("unknown", ContainerType::unknown)
        .value("nothing", ContainerType::nothing)
        .value("any", ContainerType::any);
}//function