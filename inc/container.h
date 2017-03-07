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

    virtual void print()
    {
        std::cout << "no print function defined" << std::endl;
    }//function
    
    virtual void printType();

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
    /*
    *   it's importent to coy the here since the python garbage collection will otherwise delete the data mid execution
    *   the copy call will copy the CONTAINER not the actual data (python wil then delete the old container)
    */
    void append(std::shared_ptr<Container> pC) {vElements.push_back(pC->copy());}
    bool sameTypeAs(std::shared_ptr<Container> pOther);
    void print() override
    {
        std::cout << "vector{" << std::endl;
        for(auto pElem : vElements)
        {
            if(pElem != nullptr)
                pElem->print();
            else
                std::cout << "nullptr" << std::endl;
            std::cout << "," << std::endl;
        }//for
        std::cout << "}" << std::endl;
    }//function
    
    void printType()
    {
        std::cout << "vector{";
        for(auto pElem : vElements)
        {
            if(pElem != nullptr)
                pElem->printType();
            else
                std::cout << "nullptr" << std::endl;
            std::cout << "," << std::endl;
        }//for
        std::cout << "}";
    }//function
};//class

/* function to export this module to boost python */
void exportContainer();
#endif