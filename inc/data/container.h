#include <vector>

class Container
{


};

template <Container C>
class ContainerVector<C> : public Container
{
private:
    std::vector<C> elements;

public:
    ContainerVector() 
        : 
        elements() 
    {}

    std::vector<C> elements() {return elements;}

};