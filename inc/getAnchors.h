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

    std::shared_ptr<Container> execute(std::shared_ptr<Container> pInput);

    std::shared_ptr<Container> getInputType();

    std::shared_ptr<Container> getOutputType();

};//class

void exportGetAnchors();

#endif