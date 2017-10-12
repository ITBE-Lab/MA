/** 
 * @file getAnchors.h
 * @brief Implements NlongestIntervalsAsAnchors.
 * @author Markus Schmidt
 */
#ifndef GET_ANCHORS_H
#define GET_ANCHORS_H

#include "intervalTree.h"
#include "module.h"

/**
 * @brief Extract the n longes seeds.
 * @details
 * Extracts the n longest Seed -s.
 * These can then be used as anchors for StripOfConsiderations.
 * @ingroup module
 */
class NlongestIntervalsAsAnchors : public Module
{
public:
    ///@brief number of seeds to extract
    unsigned int uiN;

    NlongestIntervalsAsAnchors(unsigned int uiN = 50)
            :  
        uiN(uiN)
    {}//constructor

    std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> vpInput);

    std::vector<ContainerType> getInputType();

    ContainerType getOutputType();

};//class

/**
 * @brief Exposes the Alignment container to boost python.
 * @ingroup export
 */
void exportGetAnchors();

#endif