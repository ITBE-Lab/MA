/** 
 * @file getAnchors.h
 * @brief Implements GetAnchors.
 * @author Markus Schmidt
 */
#ifndef GET_ANCHORS_H
#define GET_ANCHORS_H

#include "segment.h"
#include "module.h"

namespace libLAuS
{
    /**
     * @brief Extract the n longes seeds.
     * @details
     * Extracts the n longest Seed -s.
     * These can then be used as anchors for StripOfConsiderations.
     * @ingroup module
     */
    class GetAnchors : public Module
    {
    public:
        ///@brief number of seeds to extract
        unsigned int uiN;
        unsigned int uiMaxHitsPerInterval;

        GetAnchors(unsigned int uiN, unsigned int uiMaxHitsPerInterval)
                :  
            uiN(uiN),
            uiMaxHitsPerInterval(uiMaxHitsPerInterval)
        {}//constructor

        std::shared_ptr<Container> execute(ContainerVector vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - SegmentVector
         * - Pack
         * - FMIndex
         */
        ContainerVector getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - Seeds
         */
        std::shared_ptr<Container> getOutputType() const;

        std::string getName() const
        {
            return "getAnchors";
        }
    };//class
}//namespace libLAuS

/**
 * @brief Exposes the Alignment container to boost python.
 * @ingroup export
 */
void exportGetAnchors();

#endif