/** 
 * @file binarySeeding.h
 * @brief Implements a segmentation algorithm.
 * @author Markus Schmidt
 */
#ifndef BINARY_SEEDING_H
#define BINARY_SEEDING_H

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
    class BinarySeeding : public Module{
    public:
        bool bLrExtension;

        /**
         * @brief The simplified extension scheme presented in our Paper.
         * @details
         * Computes Two segments for each index as follows:
         * - extend backwards first then forwards
         * - extend forwards first then backwards
         * Starts both extensions on center.
         * Returns an interval spanning the entire covered area.
         * Segments are saved in pSegmentVector.
         */
        Interval<nucSeqIndex> lrExtension(
                nucSeqIndex center,
                std::shared_ptr<FMIndex> pFM_index,
                std::shared_ptr<NucSeq> pQuerySeq,
                std::shared_ptr<SegmentVector> pSegmentVector
            );

        /**
         * @brief The extension scheme from Li et Al.
         * @details
         * Computes all non-enclosed segments overlapping center.
         * Returns an interval spanning the entire covered area.
         * Segments are saved in pSegmentVector.
         */
        Interval<nucSeqIndex> nonEnclosedExtension(
                nucSeqIndex center,
                std::shared_ptr<FMIndex> pFM_index,
                std::shared_ptr<NucSeq> pQuerySeq,
                std::shared_ptr<SegmentVector> pSegmentVector
            );

        /*
        *    does nothing if the given interval can be found entirely on the genome.
        *    if the interval cannot be found this method splits the interval in half and repeats the step with the first half,
        *    while queuing the second half as a task in the thread pool.
        */
        void procesInterval(
                Interval<nucSeqIndex> xAreaToCover,
                std::shared_ptr<SegmentVector> pSegmentVector,
                std::shared_ptr<FMIndex> pFM_index,
                std::shared_ptr<NucSeq> pQuerySeq,
                ThreadPoolAllowingRecursiveEnqueue* pxPool
            );
        
        /*
        * functions need to be static in order to enqueue them into a threadpool
        * we need to enqueue procesInterval with an object associated
        * workaround: give the object as second (since pool will give tId as first) parameter 
        */
        static void procesIntervalStatic(
                size_t uiThreadId,
                BinarySeeding *obj,
                Interval<nucSeqIndex> xAreaToCover,
                std::shared_ptr<SegmentVector> pSegmentVector,
                std::shared_ptr<FMIndex> pFM_index,
                std::shared_ptr<NucSeq> pQuerySeq,
                ThreadPoolAllowingRecursiveEnqueue* pxPool
            )
        {
            obj->procesInterval(
                    xAreaToCover,
                    pSegmentVector,
                    pFM_index,
                    pQuerySeq,
                    pxPool
                );
        }//function

    public:
        /**
         * @brief Initialize a BinarySeeding Module
         * @details
         * if bLrExtension is True our extension scheme is used,
         * otherwise the extension scheme by Li et Al. is used.
         * Our approach is faster and computes seeds of higher quality.
         * However Li et Al.s approach will increase the overall accuracy of the alignment.
         */
        BinarySeeding(bool bLrExtension = true)
                :
            bLrExtension(bLrExtension)
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
            return "BinarySeeding";
        }
    };//class

}//namespace


/**
 * @brief exports the Segmentation @ref Module "module" to python.
 * @ingroup export
 */
void exportBinarySeeding();

#endif