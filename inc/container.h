#ifndef CONTAINER_H
#define CONTAINER_H

#include <vector>
#include <boost/python.hpp>

enum ContainerType{
    fM_index,
    nucSeq,
    packedNucSeq,
    bwtIntervalList,
    bwtInterval,
    vector,
    unknown,
};//enum

class Container
{
public:
    volatile ContainerType getType(){return ContainerType::unknown;}
    volatile bool sameTypeAs(std::shared_ptr<Container> pOther)
    {
        return pOther->getType() == getType();
    }//function

};//class

class DummyContainer : public Container
{
private:
    ContainerType iType;
public:
    DummyContainer(ContainerType iType) : iType(iType) {}//constructor
    
    ContainerType getType(){return iType;}
};//class

class ContainerVector : public Container
{
private:
    std::vector<std::shared_ptr<Container>> vElements;
public:
    ContainerVector() 
        : 
        vElements() 
    {}

    ContainerType getType(){return ContainerType::vector;}
    std::vector<std::shared_ptr<Container>> elements() {return vElements;}
    bool sameTypeAs(std::shared_ptr<Container> pOther);
};//class

/* function to export this module to boost python */
void exportContainer();
#endif