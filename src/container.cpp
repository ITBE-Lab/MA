#include "container.h"

std::vector<std::string> vContainerTypeNames =
{
    "fM_index",
    "nucSeq",
    "packedNucSeq",
    "segmentList",
    "segment",
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

bool ContainerVector::sameTypeAs(std::shared_ptr<Container> pOther)
{
    if(pOther->getType() == ContainerType::any)
        return true;
    if(pOther->getType() != getType())
        return false;
    if(vElements.size() != ((ContainerVector*)pOther.get())->vElements.size())
        return false;
    auto xItt = vElements.begin();
    auto xIttOther = ((ContainerVector*)pOther.get())->vElements.begin();
    while(xItt != vElements.end())
    {
        if( !(*xItt)->sameTypeAs(*xIttOther) )
            return false;
        xItt++;
        xIttOther++;
    }
    return true;
}//function

void exportContainer()
{
    //contianer is an abstract class and should never be initialized
	boost::python::class_<Container, std::shared_ptr<Container>>("Container", boost::python::no_init)
        .def("print_", &Container::print)
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

    //export dummy contianer
	boost::python::class_<DummyContainer, boost::python::bases<Container>, std::shared_ptr<DummyContainer>>("DummyContainer", boost::python::init<ContainerType>());

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< std::shared_ptr<DummyContainer>, std::shared_ptr<Container> >(); 
    
    //export contianer vector
	boost::python::class_<ContainerVector, boost::python::bases<Container>, std::shared_ptr<ContainerVector>>("ContainerVector")
        .def("append", &ContainerVector::append);

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< std::shared_ptr<ContainerVector>, std::shared_ptr<Container> >(); 
}//function