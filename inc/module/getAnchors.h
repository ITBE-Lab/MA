/** 
 * @file getAnchors.h
 * @brief Implements GetAnchors.
 * @author Markus Schmidt
 */
#ifndef GET_ANCHORS_H
#define GET_ANCHORS_H

#include "container/segment.h"
#include "module/stripOfConsideration.h"

namespace libMABS
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
        unsigned int uiN = 10;
        /**
         * @brief the maximal ambiguity a extracted seed shall have.
         * @details
         * Will skip any segments that would lead to more than maxAmbiguity seeds.
         */
        unsigned int uiMaxAmbiguity = 10;

        GetAnchors(unsigned int uiN, unsigned int uiMaxAmbiguity)
                :  
            uiN(uiN),
            uiMaxAmbiguity(uiMaxAmbiguity)
        {}//constructor

        GetAnchors()
        {}//constructor

        std::shared_ptr<Container> execute(std::shared_ptr<ContainerVector> vpInput);

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
}//namespace libMABS

/**
 * @brief Exposes the Alignment container to boost python.
 * @ingroup export
 */
void exportGetAnchors();

#endif