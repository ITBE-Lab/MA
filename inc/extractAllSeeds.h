/** 
 * @file extractAllSeeds.h
 * @brief Implements a Bucketing.
 * @author Markus Schmidt
 */
#ifndef EXTRACT_ALL_SEEDS_H
#define EXTRACT_ALL_SEEDS_H

#include "segmentList.h"
#include "module.h"

namespace libLAuS
{
    /**
     * @brief Used to quickly find areas with high density of @ref Seed "seeds".
     * @ingroup module
     */
    class ExtractAllSeeds: public Module
    {


    public:
        unsigned int maxNum = 10;

        ExtractAllSeeds(){}//constructor

        std::shared_ptr<Container> execute(ContainerVector vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - SegmentVector
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
            return "ExtractAllSeeds";
        }
    };//class
}//namespace libLAuS

/**
 * @brief export the ExtractAllSeeds @ref Module "module" to python.
 * @ingroup export
 */
void exportExtractAllSeeds();

#endif