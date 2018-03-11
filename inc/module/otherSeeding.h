/** 
 * @file binarySeeding.h
 * @brief Implements a segmentation algorithm.
 * @author Markus Schmidt
 */
#ifndef OTHER_SEEDING_H
#define OTHER_SEEDING_H

#include "util/system.h"
#include "module/module.h"
#include "container/segment.h"
#include "util/threadPool.h"

namespace libMA
{
    class PerfectMatch;

    /**
     * @brief Computes a maximally covering set of seeds.
     * @details
     * Can use either the extension scheme by Li et Al. or ours.
     * @ingroup module
     */
    class OtherSeeding : public Module{
    public:
        bool bBowtie;

        void bowtieExtension(
                std::shared_ptr<FMIndex> pFM_index,
                std::shared_ptr<NucSeq> pQuerySeq,
                std::shared_ptr<SegmentVector> pSegmentVector
            );

        void doBlasrExtension(
                std::shared_ptr<FMIndex> pFM_index,
                std::shared_ptr<NucSeq> pQuerySeq,
                std::shared_ptr<SegmentVector> pSegmentVector
            );

    public:
        /**
         * @brief Initialize a OtherSeeding Module
         */
        OtherSeeding(bool bBowtie = true)
                :
            bBowtie(bBowtie)
        {}//constructor
        
        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - FMIndex
         * - NucSeq
         */
        ContainerVector EXPORTED getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - SegmentVector
         */
        std::shared_ptr<Container> EXPORTED getOutputType() const;

        std::string getName() const
        {
            return "OtherSeeding";
        }
    };//class

}//namespace


/**
 * @brief exports the Segmentation @ref Module "module" to python.
 * @ingroup export
 */
void exportOtherSeeding();

#endif