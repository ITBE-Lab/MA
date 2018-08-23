/** 
 * @file mappingQuality.h
 * @brief Computes the mapping quality of Alignments.
 * @author Markus Schmidt
 */
#ifndef MAPPING_QUALITY_H
#define MAPPING_QUALITY_H

#include "module/module.h"
#include "container/alignment.h"

namespace libMA
{

    /**
     * @brief Sets the mapping quality on alignment
     * @ingroup module
     * @details
     * Given a vector of alignments this module computes the mapping quality for the
     * first alignment on the basis of the second
     * @note the name quality is missleading it rather is a mapping confidence
     */
    class MappingQuality: public Module
    {
    public:
        unsigned int uiReportNBest = defaults::uiReportN;
        float fMinSecScoreRatio = defaults::fMinSecScoreRatio;

        MappingQuality()
        {}//constructor

        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - NucSeq
         * - ContainerVector(Alignment)
         */
        ContainerVector EXPORTED getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - ContainerVector(Alignment)
         */
        std::shared_ptr<Container> EXPORTED getOutputType() const;

        std::string getName() const
        {
            return "MappingQuality";
        }

        std::string getFullDesc() const
        {
            return "MappingQuality(" + 
                std::to_string(uiReportNBest) + "," +
                std::to_string(fMinSecScoreRatio) + ")";
        }//function
    };//class
}//namspace libMA

#ifdef WITH_PYTHON
/**
 * @brief export the MappingQuality @ref Module "module" to python.
 * @ingroup export
 */
void exportMappingQuality();
#endif


#endif
