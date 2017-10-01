/** 
 * @file segmentation.h
 * @brief Implements the segmentation algorithm.
 * @author Markus Schmidt
 */
#ifndef SEGMENTATION_H
#define SEGMENTATION_H

#define DEBUG_ENABLED

#include <system.h>
#include <mutex>
#include "module.h"
#include "intervalTree.h"
#include "threadPool.h"


class PerfectMatch;


typedef DoublyLinkedList<SegmentTreeInterval>::Iterator SegTreeItt;

class Segmentation{
public:
	std::shared_ptr<SegmentTree> pSegmentTree;
private:
	std::shared_ptr<FM_Index> pxFM_index, pxRev_FM_Index;
	std::shared_ptr<NucleotideSequence> pxQuerySeq;
	bool bBreakOnAmbiguousBase;
	std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence;

	/*
	*	performs backwards extension in the given interval
	*returns the position where the backwards extension breaks
	*/
	bool canExtendFurther(std::shared_ptr<SegmentTreeInterval> pxNode, nucSeqIndex uiCurrIndex, bool bBackwards, nucSeqIndex uiQueryLength);
	/* perform forward or backwards extension (depending on bBackwards) int the given interval pxNode
	*  starts at uiStartIndex and will save any matches longer than uiMinIntervalSize in pxNode if the
	*  current extension could reach further than uiOnlyRecordHitsFurtherThan
	*/
	SaSegment extend(std::shared_ptr<SegmentTreeInterval> pxNode, nucSeqIndex uiStartIndex, bool bBackwards);
	/*
	*	does nothing if the given interval can be found entirely on the genome.
	*	if the interval cannot be found this method splits the interval in half and repeats the step with the first half,
	*	while queuing the second half as a task in the thread pool.
	*/
	void procesInterval(size_t uiThreadId, SegTreeItt pxNode, ThreadPoolAllowingRecursiveEnqueues *pxPool);


public:
	Segmentation(std::shared_ptr<FM_Index> pxFM_index, std::shared_ptr<FM_Index> pxRev_FM_index, std::shared_ptr<NucleotideSequence> pxQuerySeq,
		bool bBreakOnAmbiguousBase,
		std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence
		)
		:
		pSegmentTree(new SegmentTree(pxQuerySeq->length())),
		pxFM_index(pxFM_index),
		pxRev_FM_Index(pxRev_FM_index),
		pxQuerySeq(pxQuerySeq),
		bBreakOnAmbiguousBase(bBreakOnAmbiguousBase),
		pxRefSequence(pxRefSequence)
	{}//constructor

	void segment();
};//class

class SegmentationContainer : public Module
{
public:
	bool bBreakOnAmbiguousBase;
	unsigned int uiMaxHitsPerInterval;

	SegmentationContainer
	(
		bool bBreakOnAmbiguousBase = true,
		unsigned int uiMaxHitsPerInterval = 10000
	)
			:
		bBreakOnAmbiguousBase(bBreakOnAmbiguousBase),
		uiMaxHitsPerInterval(uiMaxHitsPerInterval)
	{}//constructor

	std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> vpInput);

	std::vector<ContainerType> getInputType();
    std::vector<ContainerType> getOutputType();
};//class

void exportSegmentation();

#endif