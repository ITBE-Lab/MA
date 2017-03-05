#ifndef CONTAINER_H
#define CONTAINER_H

#include <vector>
#include <boost/python.hpp>

class Container
{


};

template<class C>
class ContainerVector : public Container
{
private:
    std::vector<C> vElements;

public:
    ContainerVector() 
        : 
        vElements() 
    {}

    std::vector<C> elements() {return vElements;}

};

/* function to export this module to boost python */
void exportContainer();
#endif