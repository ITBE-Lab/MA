#ifndef INTERVALTREE_H
#define INTERVALTREE_H

#include <memory>
#include <assert.h>
#include "string"
#include <iostream>
#include <mutex> 
#include <vector>
#include "FM_index.h"
#include <list>
#include <thread>
#include <boost/python.hpp>
#include "doublyLinkedList.h"
#include "seed.h"

#define confMETA_MEASURE_DURATION ( 1 )



class SaSegment: public Container, public Interval<nucSeqIndex> {
private:
	SA_IndexInterval xSaInterval;
	bool bForw;
public:
	/**
	*	start and end on querry
	*	sa index interval holds the bwt interval
	*/
	SaSegment(nucSeqIndex uiStart, nucSeqIndex uiSize, SA_IndexInterval xSaInterval, bool bForw)
			:
		Interval(uiStart, uiSize),
		xSaInterval(xSaInterval),
		bForw(bForw)
	{}//constructor

	ContainerType getType() const {return ContainerType::segment;}
	const SA_IndexInterval& saInterval() const
	{
		return xSaInterval;
	}//function
	bool isForward() const
	{
		return bForw;
	}//function
}; // class ( Segment )

class SegmentTreeInterval: public Container, public Interval<nucSeqIndex>
{
private:
	/** list of the perfect matches found through backwards / forward extension */
	std::list<SaSegment> lxSaSegment;
	/** list of the longest perfect matches found through backwards / forward extension */
	std::list<SaSegment> lxSaAnchorSegment;

public:
	SegmentTreeInterval(const nucSeqIndex uiStart, const nucSeqIndex uiSize)
		:
		Interval(uiStart, uiSize),
		lxSaSegment(),
		lxSaAnchorSegment()
	{}//constructor
	
	ContainerType getType(){return ContainerType::segment;}//function


	/* prints information about this node; thread save */
	void print(std::ostream& xOs) const
	{
		xOs << "(" << std::to_string(this->start()) << "," << std::to_string(this->end()) << ")";
	}//function
	/* push back an interval of perfect matches
	 * the interval contains uiLengthInBwt individual perfect matches of (uiStartOfIntervalOnQuery, uiEndOfIntervalOnQuery) on the reference sequence
	 * bForwHit is required because we need to extract the starting positions of every single match in the interval using the same fm_index used to calculate the interval
	 * since we use 2 fm_indecies one for forward one for reverse we need to remember which of the two was used
	*/
	void push_back(SaSegment interval, bool bAnchor);

	nucSeqIndex getCenter() const 
	{
		return start() + size() / 2; 
	}//function

	nucSeqIndex getValue() const 
	{
		return size(); 
	}//function

	/* calls fDo for all recorded hits.
	*  Note that pushBackBwtInterval records an interval of hits
	*  ulHit: the position of the hit on the reference sequence
	*  uiQueryBegin: the position of the hit on the query sequence
	*  uiQueryEnd: the end of the hit on the query sequence
	*  bForwHit: true if this hit was produced through forward extension
	*/
	void forEachSeed(
			std::shared_ptr<FM_Index> pxFM_Index,
			std::shared_ptr<FM_Index> pxRev_FM_Index,
			unsigned int uiMaxNumHitsPerInterval,
			bool bSkipLongerIntervals,
			bool bAnchorOnly,
			std::function<void(Seed s)> fDo
		)
	{
		//iterate over all the intervals that have been recorded using pushBackBwtInterval()
		for (SaSegment xSegment : bAnchorOnly ? lxSaAnchorSegment : lxSaSegment)
		{
			//if the interval contains more than uiMaxNumHitsPerInterval hits it's of no importance and will produce nothing but noise

			//if bSkipLongerIntervals is not set uiJump by is used to not return more than 
			//uiMaxNumHitsPerInterval
			t_bwtIndex uiJumpBy = 1;
			if (xSegment.saInterval().size() > uiMaxNumHitsPerInterval)
			{
				if (bSkipLongerIntervals)
					continue;
				uiJumpBy = xSegment.saInterval().size() / uiMaxNumHitsPerInterval; 
			}//if

			//if the hit was generated using the reversed fm_index we should use the according fm_index in order to 
			//extract the index of the hit on the reference sequence. same for the forward fm_index
			std::shared_ptr<FM_Index> pxUsedFmIndex;
			if (xSegment.isForward())
				pxUsedFmIndex = pxRev_FM_Index;
			else
				pxUsedFmIndex = pxFM_Index;

			//iterate over the interval in the BWT
			for (
					auto ulCurrPos = xSegment.saInterval().start(); 
					ulCurrPos < xSegment.saInterval().end(); 
					ulCurrPos += uiJumpBy
				)
			{
				//calculate the referenceIndex using pxUsedFmIndex->bwt_sa() and call fDo for every match individually
				auto ulIndexOnRefSeq = pxUsedFmIndex->bwt_sa(ulCurrPos);
				/* if the match was calculated using the fm-index of the reversed sequence:
				 * we acquire the index of the beginning of the match on the reversed sequence by calling bwt_sa()
				 * but we actually want the beginning of the match on the normal sequence, so we need to subtract the END of the match from the reference sequence length
				 */
				if (xSegment.isForward())
					ulIndexOnRefSeq = pxUsedFmIndex->getRefSeqLength() - (ulIndexOnRefSeq + xSegment.size());
				//call the given function
				fDo(Seed(xSegment.start(), xSegment.size(), ulIndexOnRefSeq));
			}//for
		}//for
	}//function

	std::vector<std::shared_ptr<NucleotideSequence>> getRefHits(
			std::shared_ptr<FM_Index> pxFM_Index, 
			std::shared_ptr<FM_Index> pxRev_FM_Index,
			std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefPack
		)
	{
		std::vector<std::shared_ptr<NucleotideSequence>> vpRet = 
			std::vector<std::shared_ptr<NucleotideSequence>>();
		forEachSeed(
			pxFM_Index,
			pxRev_FM_Index,
			100000,
			false,
			false,
			[&](Seed xS)
			{
				vpRet.push_back(pxRefPack->vExtract(xS.start_ref(), xS.end_ref()));
			}//lambda
		);//forall
		return vpRet;
	}//function
};//class

class SegmentTree : public DoublyLinkedList<SegmentTreeInterval>, public Container{

public:
	/*
	*	sets up the interval tree with two leaves and one initial interval comprising the whole query
	*	note that the tree is internally represented as a DoublyLinkedList since only the leaves are of relevance
	*/
	SegmentTree(const nucSeqIndex uiQuerryLength)
		:
		DoublyLinkedList()
	{
		std::shared_ptr<SegmentTreeInterval> pxRoot(new SegmentTreeInterval(0, uiQuerryLength));
		push_back(pxRoot);
	}//constructor
	SegmentTree()
	{}//constructor
	
	ContainerType getType(){return ContainerType::segmentList;}//function

	/*not thread save; prints basic information about the segment tree*/
	void print(std::ostream &xOut) const
	{
		forEach([&xOut](std::shared_ptr<SegmentTreeInterval> pxNode){ pxNode->print(xOut); });
	}//function
};

std::ostream& operator<<(std::ostream& xOs, const SegmentTree& rxTree);
std::ostream& operator<<(std::ostream& xOs, const SegmentTreeInterval &rxNode);


void exportIntervalTree();


#endif