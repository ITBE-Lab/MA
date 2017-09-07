#ifndef SEGMENTATION_H
#define SEGMENTATION_H


#include "intervalTree.h"
#include "FM_index.h"
#include <condition_variable>
#include "threadPool.h"
#include "graphicalMethod.h"
#include "module.h"
#include "container.h"
#include <system.h>


class PerfectMatch;

typedef DoublyLinkedList<SegmentTreeInterval>::Iterator SegTreeItt;

class Segmentation{
public:
	std::shared_ptr<SegmentTree> pSegmentTree;
private:
	std::shared_ptr<FM_Index> pxFM_index, pxRev_FM_Index;
	std::shared_ptr<NucleotideSequence> pxQuerySeq;
	//std::shared_ptr<AnchorMatchList> pxAnchorMatchList;  TODO: move me into my own module
	bool bBreakOnAmbiguousBase;
	bool bSkipLongBWTIntervals;
	nucSeqIndex uiMinIntervalSize;
	std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence;
	unsigned int uiMaxHitsPerInterval;
#if confGENEREATE_ALIGNMENT_QUALITY_OUTPUT
	AlignmentQuality *pxQuality;
#endif

	/*
	*	performs backwards extension in the given interval
	*returns the position where the backwards extension breaks
	*/
	bool canExtendFurther(std::shared_ptr<SegmentTreeInterval> pxNode, nucSeqIndex uiCurrIndex, bool bBackwards, nucSeqIndex uiQueryLength);
	/* perform forward or backwards extension (depending on bBackwards) int the given interval pxNode
	*  starts at uiStartIndex and will save any matches longer than uiMinIntervalSize in pxNode if the
	*  current extension could reach further than uiOnlyRecordHitsFurtherThan
	*/
	nucSeqIndex extend(std::shared_ptr<SegmentTreeInterval> pxNode, nucSeqIndex uiStartIndex, bool bBackwards, nucSeqIndex uiOnlyRecordHitsFurtherThan);
	/*
	*	does nothing if the given interval can be found entirely on the genome.
	*	if the interval cannot be found this method splits the interval in half and repeats the step with the first half,
	*	while queuing the second half as a task in the thread pool.
	*/
	void procesInterval(size_t uiThreadId, SegTreeItt pxNode, ThreadPoolAllowingRecursiveEnqueues *pxPool);

	/*
	*	deprecated since the anchor matches will be extracted after the segmentation process.
	* 	the finding achor matches process will have it's own module
	*/
	//void saveHits(std::shared_ptr<SegmentTreeInterval> pxNode, size_t uiThreadId);
	void forEachNonBridgingHitOnTheRefSeq(std::shared_ptr<SegmentTreeInterval> pxNode, bool bAnchorOnly,
		std::function<void(nucSeqIndex ulIndexOnRefSeq, nucSeqIndex uiQueryBegin, nucSeqIndex uiQueryEnd)> fDo);
	void forEachNonBridgingPerfectMatch(std::shared_ptr<SegmentTreeInterval> pxNode, bool bAnchorOnly,
		std::function<void(std::shared_ptr<PerfectMatch>)> fDo);

public:
	Segmentation(std::shared_ptr<FM_Index> pxFM_index, std::shared_ptr<FM_Index> pxRev_FM_index, std::shared_ptr<NucleotideSequence> pxQuerySeq,
		bool bBreakOnAmbiguousBase, bool bSkipLongBWTIntervals, unsigned int uiMinIntervalSize, unsigned int uiMaxHitsPerInterval,
		std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence
#if confGENEREATE_ALIGNMENT_QUALITY_OUTPUT
		,AlignmentQuality *pxQuality
#endif
		)
		:
		pSegmentTree(new SegmentTree(pxQuerySeq->length() - 2)),
		pxFM_index(pxFM_index),
		pxRev_FM_Index(pxRev_FM_index),
		pxQuerySeq(pxQuerySeq),
		bBreakOnAmbiguousBase(bBreakOnAmbiguousBase),
		bSkipLongBWTIntervals(bSkipLongBWTIntervals),
		uiMinIntervalSize(uiMinIntervalSize),
		pxRefSequence(pxRefSequence),
		uiMaxHitsPerInterval(uiMaxHitsPerInterval)
#if confGENEREATE_ALIGNMENT_QUALITY_OUTPUT
		,pxQuality(pxQuality)
#endif
	{}//constructor

	void segment();
};//class

class SegmentationContainer : public Module
{
public:
	bool bBreakOnAmbiguousBase;
	bool bSkipLongBWTIntervals;
	nucSeqIndex uiMinIntervalSize;
	unsigned int uiMaxHitsPerInterval;

	SegmentationContainer
	(
		bool bBreakOnAmbiguousBase = true,
		bool bSkipLongBWTIntervals = true,
		nucSeqIndex uiMinIntervalSize = 10,
		unsigned int uiMaxHitsPerInterval = 10000
	)
			:
		bBreakOnAmbiguousBase(bBreakOnAmbiguousBase),
		bSkipLongBWTIntervals(bSkipLongBWTIntervals),
		uiMinIntervalSize(uiMinIntervalSize),
		uiMaxHitsPerInterval(uiMaxHitsPerInterval)
	{}//constructor

	std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> vpInput);

	std::vector<ContainerType> getInputType();
    std::vector<ContainerType> getOutputType();
};//class

void exportSegmentation();

#endif