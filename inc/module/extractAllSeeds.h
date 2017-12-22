/** 
 * @file extractAllSeeds.h
 * @brief Implements a Bucketing.
 * @author Markus Schmidt
 */
#ifndef EXTRACT_ALL_SEEDS_H
#define EXTRACT_ALL_SEEDS_H

#include "container/segment.h"
#include "module/module.h"

namespace libMABS
{
    /**
     * @brief Extracts all Seeds from a SegmentList.
     * @ingroup module
     */
    class ExtractAllSeeds: public Module
    {
    public:
        /**
         * @brief the maximal ambiguity a extracted seed shall have.
         * @details
         * Will skip any segments that would lead to more than maxAmbiguity seeds.
         */
        unsigned int maxAmbiguity = 10;

        ExtractAllSeeds(){}//default constructor

        ExtractAllSeeds(unsigned int maxAmbiguity)
                :
            maxAmbiguity(maxAmbiguity)
        {}//constructor

        std::shared_ptr<Container> execute(std::shared_ptr<ContainerVector> vpInput);

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
}//namespace libMABS

/**
 * @brief export the ExtractAllSeeds @ref Module "module" to python.
 * @ingroup export
 */
void exportExtractAllSeeds();

#endif