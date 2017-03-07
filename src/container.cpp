#include "container.h"

bool ContainerVector::sameTypeAs(std::shared_ptr<Container> pOther)
{
    if(pOther->getType() != getType())
        return false;
    if(elements().size() != ((ContainerVector*)pOther.get())->elements().size())
        return false;
    auto xItt = elements().begin();
    auto xIttOther = ((ContainerVector*)pOther.get())->elements().begin();
    while(xItt != elements().end())
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
	boost::python::class_<Container, std::shared_ptr<Container>>("Container", boost::python::no_init);
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