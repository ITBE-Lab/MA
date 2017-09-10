#ifndef GET_ANCHORS_H
#define GET_ANCHORS_H

#include "boost/python.hpp"
#include "intervalTree.h"
#include "module.h"

class NlongestIntervalsAsAnchors : public Module
{
public:
    unsigned int uiN;

    NlongestIntervalsAsAnchors(unsigned int uiN = 50)
            :  
        uiN(uiN)
    {}//constructor

    std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> vpInput);

    std::vector<ContainerType> getInputType();

    std::vector<ContainerType> getOutputType();

};//class

void exportGetAnchors();

#endif