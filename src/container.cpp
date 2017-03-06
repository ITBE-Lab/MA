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
} 

void exportContainer()
{
    //contianer is an abstract class and should never be initialized
	boost::python::class_<Container, std::shared_ptr<Container>>("Container", boost::python::no_init);

    //export the containertype enum
    boost::python::enum_<ContainerType>("ContainerType")
        .value("unknown", ContainerType::unknown)
        .value("vector", ContainerType::vector)
        .value("FM_index", ContainerType::fM_index)
        .value("nucSeq", ContainerType::nucSeq)
        .value("bwtIntervalList", ContainerType::bwtIntervalList)
        .value("bwtInterval", ContainerType::bwtInterval)
        .value("packedNucSeq", ContainerType::packedNucSeq);

    //export dummy contianer
	boost::python::class_<DummyContainer, boost::python::bases<Container>, std::shared_ptr<DummyContainer>>("DummyContainer", boost::python::init<ContainerType>());
    
	//tell boost python that it's possible to convert shared pointers with these classes
    boost::python::implicitly_convertible<std::shared_ptr<DummyContainer>,std::shared_ptr<Container>>();
}