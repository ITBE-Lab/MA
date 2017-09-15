#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

#include "boost/python.hpp"
#include "seed.h"
#include "alignment.h"
#include "graphicalMethod.h"
#include <memory>
#include <iostream>
#include "debug.h"

class NeedlemanWunsch : public Module
{
public:

    NeedlemanWunsch()
    {}//constructor

    std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> vpInput);

    std::vector<ContainerType> getInputType();

    std::vector<ContainerType> getOutputType();

};//class

void exportNeedlemanWunsch();

#endif