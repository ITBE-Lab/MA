/** 
 * @file getAnchors.h
 * @brief Implements NlongestIntervalsAsAnchors.
 * @author Markus Schmidt
 */
#ifndef GET_ANCHORS_H
#define GET_ANCHORS_H

#include "segmentList.h"
#include "cppModule.h"

/**
 * @brief Extract the n longes seeds.
 * @details
 * Extracts the n longest Seed -s.
 * These can then be used as anchors for StripOfConsiderations.
 * @ingroup module
 */
class NlongestIntervalsAsAnchors : public CppModule
{
public:
    ///@brief number of seeds to extract
    unsigned int uiN;
    unsigned int uiMaxHitsPerInterval;

    NlongestIntervalsAsAnchors(unsigned int uiN, unsigned int uiMaxHitsPerInterval)
            :  
        uiN(uiN),
        uiMaxHitsPerInterval(uiMaxHitsPerInterval)
    {}//constructor

    std::shared_ptr<Container> execute(ContainerVector vpInput);

    ContainerVector getInputType() const;

    std::shared_ptr<Container> getOutputType() const;

    std::string getName() const
    {
        return "getAnchors";
    }
};//class

/**
 * @brief Exposes the Alignment container to boost python.
 * @ingroup export
 */
void exportGetAnchors();

#endif