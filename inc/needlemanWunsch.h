#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

#define DEBUG_LEVEL 1

#include "boost/python.hpp"
#include "intervalTree.h"
#include "module.h"
#include "sequence.h"
#include "aligner.h"
#include <memory>
#include <iostream>
#include "meta_programming.h"

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