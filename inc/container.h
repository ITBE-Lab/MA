#ifndef CONTAINER_H
#define CONTAINER_H

#include <vector>
#include <boost/python.hpp>
#include <iostream>

enum ContainerType{
    fM_index,
    nucSeq,
    packedNucSeq,
    segmentList,
    segment,
    stripOfConsideration,
    stripOfConsiderationList,
    vector,
    unknown,
    nothing,
    any,
};//enum


class Container
{
public:
    virtual ContainerType getType(){return ContainerType::unknown;}

    virtual bool sameTypeAs(std::shared_ptr<Container> pOther)
    {
        if(pOther->getType() == ContainerType::any)
            return true;
        if(getType() == ContainerType::any)
            return true;
        return pOther->getType() == getType();
    }//function

    virtual std::string toString()
    {
        return std::string("no print function defined");
    }//function
    
    virtual std::string getTypeInfo();

    //this should always be overwritten
    virtual std::shared_ptr<Container> copy()
    {
        return std::shared_ptr<Container>();
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

/* function to export this module to boost python */
void exportContainer();
#endif