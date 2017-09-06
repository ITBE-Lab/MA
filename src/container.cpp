#include "container.h"

std::vector<std::string> vContainerTypeNames =
{
    "fM_index",
    "nucSeq",
    "packedNucSeq",
    "segmentList",
    "segment",
    "stripOfConsideration",
    "stripOfConsiderationList",
    "vector",
    "unknown",
    "nothing",
    "any",
};//array

std::string Container::getTypeInfo()
{
    if(vContainerTypeNames.size() > getType())
        return std::string(vContainerTypeNames[getType()]);
    return std::string("unknown");
}//function

void exportContainer()
{
    //contianer is an abstract class and should never be initialized
	boost::python::class_<Container, std::shared_ptr<Container>>("Container", boost::python::no_init)
        .def("getTypeInfo", &Container::getTypeInfo);
    boost::python::register_ptr_to_python< std::shared_ptr<Container> >();


    //export the containertype enum
    boost::python::enum_<ContainerType>("ContainerType")
        .value("unknown", ContainerType::unknown)
        .value("vector", ContainerType::vector)
        .value("FM_index", ContainerType::fM_index)
        .value("nucSeq", ContainerType::nucSeq)
        .value("segmentList", ContainerType::segmentList)
        .value("segment", ContainerType::segment)
        .value("packedNucSeq", ContainerType::packedNucSeq)
        .value("nothing", ContainerType::nothing);
}//function