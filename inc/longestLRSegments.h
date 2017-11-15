/** 
 * @file longestNonEnclosedSegments.h
 * @brief Implements a segmentation algorithm.
 * @author Markus Schmidt
 */
#ifndef LONGEST_L_R_SEGMENTS_H
#define LONGEST_L_R_SEGMENTS_H

#define DEBUG_ENABLED

#include "system.h"
#include "cppModule.h"
#include "intervalTree.h"
#include "threadPool.h"


class PerfectMatch;

/**
 * @brief Computes a set of maximal non-enclosed seeds.
 * @ingroup module
 */
class LongestLRSegments : public CppModule{
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
	static SaSegment extend(
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
			DoublyLinkedList<SegmentTreeInterval>::Iterator pxNode, 
			std::shared_ptr<SegmentTree> pSegmentTree,
			std::shared_ptr<FM_Index> pFM_index,
			std::shared_ptr<NucleotideSequence> pQuerySeq,
			ThreadPoolAllowingRecursiveEnqueues* pxPool
		);


public:
	
	std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> vpInput);

	std::vector<std::shared_ptr<Container>> getInputType() const;
    std::shared_ptr<Container> getOutputType() const;

    std::string getName() const
    {
        return "LongestLRSegments";
    }
};//class


/**
 * @brief exports the Segmentation @ref CppModule "module" to python.
 * @ingroup export
 */
void exportLongestLRSegments();

#endif