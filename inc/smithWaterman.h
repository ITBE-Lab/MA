#ifndef SMITH_WATERMAN_H
#define SMITH_WATERMAN_H

#include "boost/python.hpp"
#include "intervalTree.h"
#include "module.h"

class SmithWaterman : public Module
{
public:

    SmithWaterman()
    {}//constructor

    std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> vpInput);

    std::vector<ContainerType> getInputType();

    std::vector<ContainerType> getOutputType();

};//class

void exportSmithWaterman();

#endif