#ifndef CONTAINER_H
#define CONTAINER_H

#include <vector>

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


#endif