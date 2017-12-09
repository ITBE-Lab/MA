/** 
 * @file longestNonEnclosedSegments.h
 * @brief Implements a segmentation algorithm.
 * @author Markus Schmidt
 */
#ifndef LONGEST_NON_ENCLOSED_SEGMENTS_H
#define LONGEST_NON_ENCLOSED_SEGMENTS_H

#define DEBUG_ENABLED

#include "system.h"
#include "cppModule.h"
#include "intervalTree.h"
#include "threadPool.h"


class PerfectMatch;

/**
 * @brief Computes all non-enclosed seeds.
 * @details
 * This is the the BWA seeding strategy, parallelized using binary seeding.
 * Gives you all non-enclosed seeds.
 * Therefore there are a lot of bad quality seeds and strict filters need to be applied.
 * @ingroup module
 */
class LongestNonEnclosedSegments : public CppModule{
private:

	static SA_IndexInterval extend_backward(
			const SA_IndexInterval &ik, 
			const uint8_t c, 
			std::shared_ptr<FM_Index> pFM_index
		);

	/* perform forward or backwards extension (depending on bBackwards) int the given interval pxNode
	*  starts at uiStartIndex and will save any matches longer than uiMinIntervalSize in pxNode if the
	*  current extension could reach further than uiOnlyRecordHitsFurtherThan
	*/
	static Interval<nucSeqIndex> extend(
			std::shared_ptr<SegmentTreeInterval> pxNode,
			std::shared_ptr<FM_Index> pFM_index,
			std::shared_ptr<NucleotideSequence> pQuerySeq
		);
	/*
	*	does nothing if the given interval can be found entirely on the genome.
	*	if the interval cannot be found this method splits the interval in half and repeats the step with the first half,
	*	while queuing the second half as a task in the thread pool.
	*/
	static void procesInterval(
			size_t uiThreadId, 
			SegmentTree::iterator pxNode, 
			std::shared_ptr<SegmentTree> pSegmentTree,
			std::shared_ptr<FM_Index> pFM_index,
			std::shared_ptr<NucleotideSequence> pQuerySeq,
			ThreadPoolAllowingRecursiveEnqueues* pxPool
		);


public:
	
	std::shared_ptr<Container> execute(ContainerVector vpInput);

	ContainerVector getInputType() const;
    std::shared_ptr<Container> getOutputType() const;

    std::string getName() const
    {
        return "LongestNonEnclosedSegments";
    }
};//class


/**
 * @brief exports the Segmentation @ref CppModule "module" to python.
 * @ingroup export
 */
void exportLongestNonEnclosedSegments();

#endif