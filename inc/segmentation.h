#ifndef SEGMENTATION_H
#define SEGMENTATION_H


#include "intervalTree.h"
#include "FM_index.h"
#include <condition_variable>
#include "GraphicalMethod.h"
#include <system.h>



typedef DoublyLinkedList<SegmentTreeInterval>::Iterator SegTreeItt;

class Segmentation{
private:
	SegmentTree xSegmentTree;
	std::shared_ptr<FM_Index> pxFM_index, pxRev_FM_Index;
	std::shared_ptr<NucleotideSequence> pxQuerySeq;
	nucSeqIndex uiMinIntervalSize;
	bool bBreakOnAmbiguousBase;
	bool bSkipLongBWTIntervals;
	std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence;
	unsigned int uiMaxHitsPerInterval;
	unsigned int uiNumSegmentsAsAnchors;
	std::shared_ptr<AnchorMatchList> pxAnchorMatchList;
#if confGENEREATE_ALIGNMENT_QUALITY_OUTPUT
	AlignmentQuality *pxQuality;
#endif

	/*
	*	performs backwards extension in the given interval
	*	returns the position where the backwards extension breaks
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


	void saveHits(std::shared_ptr<SegmentTreeInterval> pxNode, size_t uiThreadId);
	void forEachNonBridgingHitOnTheRefSeq(std::shared_ptr<SegmentTreeInterval> pxNode, bool bAnchorOnly,
		std::function<void(nucSeqIndex ulIndexOnRefSeq, nucSeqIndex uiQueryBegin, nucSeqIndex uiQueryEnd)> fDo);
	void forEachNonBridgingPerfectMatch(std::shared_ptr<SegmentTreeInterval> pxNode, bool bAnchorOnly,
		std::function<void(std::shared_ptr<PerfectMatch>)> fDo);

public:
	Segmentation(std::shared_ptr<FM_Index> pxFM_index, std::shared_ptr<FM_Index> pxRev_FM_index, std::shared_ptr<NucleotideSequence> pxQuerySeq,
		bool bBreakOnAmbiguousBase, bool bSkipLongBWTIntervals, unsigned int uiMinIntervalSize, unsigned int uiMaxHitsPerInterval, unsigned int uiNumSegmentsAsAnchors,
		std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence, std::shared_ptr<AnchorMatchList> pxAnchorMatchList
#if confGENEREATE_ALIGNMENT_QUALITY_OUTPUT
		,AlignmentQuality *pxQuality
#endif
		)
		:
		pxFM_index(pxFM_index),
		pxRev_FM_Index(pxRev_FM_index),
		xSegmentTree(pxQuerySeq->length() - 2),
		pxQuerySeq(pxQuerySeq),
		bBreakOnAmbiguousBase(bBreakOnAmbiguousBase),
		bSkipLongBWTIntervals(bSkipLongBWTIntervals),
		uiMinIntervalSize(uiMinIntervalSize),
		pxRefSequence(pxRefSequence),
		uiMaxHitsPerInterval(uiMaxHitsPerInterval),
		uiNumSegmentsAsAnchors(uiNumSegmentsAsAnchors),
		pxAnchorMatchList(pxAnchorMatchList)
#if confGENEREATE_ALIGNMENT_QUALITY_OUTPUT
		,pxQuality(pxQuality)
#endif
	{}//constructor

	/*
	*	performs the entire segmentation process
	*/ 
	void segment();

	void printSegmentTree(std::ostream &xOut) { xSegmentTree.print(xOut); xOut << std::endl; }//function

};


#endif