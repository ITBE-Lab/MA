#ifndef CONTAINER_H
#define CONTAINER_H

#include <vector>
#include <boost/python.hpp>

enum ContainerType{
    FM_index,
    nucSeq,
    vector,
    unknown,
};

class Container
{
public:
    volatile ContainerType getType(){return ContainerType::unknown;}
    volatile bool sameTypeAs(std::shared_ptr<Container> &xOther)
    {
        return xOther->getType() == getType();
    }

};

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
    bool sameTypeAs(std::shared_ptr<Container> pOther)
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
};

/* function to export this module to boost python */
void exportContainer();
#endif