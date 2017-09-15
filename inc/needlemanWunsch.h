/** 
 * @file needlemanWunsch.h
 * @brief Implements NMW
 * @author Markus Schmidt
 */
#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

#include "alignment.h"
#include "graphicalMethod.h"
#include "module.h"

/**
 * @brief implements NMW
 * @details
 * Returns a finished alignment if given a sound selection of seeds.
 */
class NeedlemanWunsch : public Module
{
public:

    NeedlemanWunsch()
    {}//constructor

    //overload
    std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> vpInput);

    //overload
    std::vector<ContainerType> getInputType();

    //overload
    std::vector<ContainerType> getOutputType();

};//class

/**
 * @brief Exposes the Alignment container to boost python.
 */
void exportNeedlemanWunsch();

#endif