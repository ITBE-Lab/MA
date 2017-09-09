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
//#include <csignal>

//#define DEBUG_ALIGNER 
//#define DEBUG_CHECK_INTERVALS 
//changing this to 1 might break a lot of stuff (old code) [i don't think it's worth it to go through and repair it though]
#define USE_BUCKET_CLUTERING ( 0 )

#define confMETA_MEASURE_DURATION ( 1 )
#define confGENEREATE_ALIGNMENT_QUALITY_OUTPUT ( 0 )

#ifdef DEBUG_CHECK_INTERVALS
static std::mutex xMutexDebugCheckIntervals;
#endif



//any index on the query or reference nucleotide sequence is given in this datatype
typedef uint64_t nucSeqIndex;


#if confGENEREATE_ALIGNMENT_QUALITY_OUTPUT
class AlignmentQuality{
public:

	nucSeqIndex uiSeedCoverage;
	nucSeqIndex uiQueryLength;
	nucSeqIndex uiLongestGapInAlignment;
	nucSeqIndex uiLongestSeedInAlignment;
	unsigned int uiAmountOfGapsInAligment;
	unsigned int uiAmountOfSeedsInAligment;
	nucSeqIndex uiAmountOfOverlapsInAligment;

	unsigned int uiAmountSegments;
	nucSeqIndex uiLongestSegment;
	nucSeqIndex uiLongestGapBetweenSegments;
	unsigned int uiAmountAnchorSegments;
	//unused because of uiMaxHitsPerInterval
	unsigned int uiAmountOfSkippedBWAIntervals = 0;
	unsigned int uiMinIntervalSize;
	unsigned int uiNumSegmentsAsAnchors;

	double fBWADur;
	double fTotal;
	double fSegmentation;
	double fAnchorMatches;
	double fLineSweep;
	unsigned int uiBWACoverage;

	AlignmentQuality()
		:
		uiSeedCoverage(0),
		uiQueryLength(0),
		uiLongestGapInAlignment(0),
		uiLongestSeedInAlignment(0),
		uiAmountOfSeedsInAligment(0),
		uiAmountOfGapsInAligment(0),
		uiAmountOfOverlapsInAligment(0),
		uiAmountSegments(0),
		uiLongestSegment(0),
		uiLongestGapBetweenSegments(0),
		uiAmountAnchorSegments(0),
		uiAmountOfSkippedBWAIntervals(0),
		uiMinIntervalSize(0),
		uiNumSegmentsAsAnchors(0),
		fBWADur(0),
		fTotal(0),
		fSegmentation(0),
		fAnchorMatches(0),
		fLineSweep(0),
		uiBWACoverage(0)
	{}//constructor	


	static void printQualityHeader(std::ostream &xOut)
	{
		xOut << "seedCoverage \t uiQueryLength \t uiLongestGapInAlignment \t uiLongestSeedInAlignment \t"
		 	 << " uiAmountOfGapsInAligment \t uiAmountOfSeedsInAligment \t uiAmountOfOverlapsInAligment \t"
			 << " uiAmountSegments \t uiLongestSegment \t uiLongestGapBetweenSegments \t uiAmountAnchorSegments \t"
			 << " uiAmountOfSkippedBWAIntervals \t uiMinIntervalSize \t uiNumSegmentsAsAnchors \t fBWADur \t"
			 << " fTotal \t fSegmentation \t fAnchorMatches \t fLineSweep \t uiBWACoverage" << std::endl;
	}

	void printQuality(std::ostream &xOut)
	{
		xOut << uiSeedCoverage << "//seedCoverage" << std::endl;
		xOut << uiQueryLength << "//uiQueryLength" << std::endl;
		xOut << uiLongestGapInAlignment << "//uiLongestGapInAlignment" << std::endl;
		xOut << uiLongestSeedInAlignment << "//uiLongestSeedInAlignment" << std::endl;
		xOut << uiAmountOfGapsInAligment << "//uiAmountOfGapsInAligment" << std::endl;
		xOut << uiAmountOfSeedsInAligment << "//uiAmountOfSeedsInAligment" << std::endl;
		xOut << uiAmountOfOverlapsInAligment << "//uiAmountOfOverlapsInAligment" << std::endl;
		xOut << uiAmountSegments << "//uiAmountSegments" << std::endl;
		xOut << uiLongestSegment << "//uiLongestSegment" << std::endl;
		xOut << uiLongestGapBetweenSegments << "//uiLongestGapBetweenSegments" << std::endl;
		xOut << uiAmountAnchorSegments << "//uiAmountAnchorSegments" << std::endl;
		xOut << uiAmountOfSkippedBWAIntervals << "//uiAmountOfSkippedBWAIntervals" << std::endl;
		xOut << uiMinIntervalSize << "//uiMinIntervalSize" << std::endl;
		xOut << uiNumSegmentsAsAnchors << "//uiNumSegmentsAsAnchors" << std::endl;
		xOut << fBWADur << "//fBWADur" << std::endl;
		xOut << fTotal << "//fTotal" << std::endl;
		xOut << fSegmentation << "//fSegmentation" << std::endl;
		xOut << fAnchorMatches << "//fAnchorMatches" << std::endl;
		xOut << fLineSweep << "//fLineSweep" << std::endl;
		xOut << uiBWACoverage << "//uiBWACoverage" << std::endl;
	}//function

	void printQualityExcelFormat(std::ostream &xOut)
	{
		xOut << uiSeedCoverage << '\t';
		xOut << uiQueryLength << '\t';
		xOut << uiLongestGapInAlignment << '\t';
		xOut << uiLongestSeedInAlignment << '\t';
		xOut << uiAmountOfGapsInAligment << '\t';
		xOut << uiAmountOfSeedsInAligment << '\t';
		xOut << uiAmountOfOverlapsInAligment << '\t';
		xOut << uiAmountSegments << '\t';
		xOut << uiLongestSegment << '\t';
		xOut << uiLongestGapBetweenSegments << '\t';
		xOut << uiAmountAnchorSegments << '\t';
		xOut << uiAmountOfSkippedBWAIntervals << '\t';
		xOut << uiMinIntervalSize << '\t'; 
		xOut << uiNumSegmentsAsAnchors << '\t';
		xOut << fBWADur << '\t';
		xOut << fTotal << '\t';
		xOut << fSegmentation << '\t';
		xOut << fAnchorMatches << '\t';
		xOut << fLineSweep << '\t';
		xOut << uiBWACoverage << std::endl;
	}//function
};
#endif



#ifdef DEBUG_INTERVALTREE
	/* when debugging each node gets its unique id in order to track it */
	static unsigned int iNextId = 0;
#endif // DEBUG_INTERVALTREE


class SegmentTreeInterval: public Container
{
private:
	struct PerfectMatch{
		t_bwtIndex uiBwtIntervalIndex;
		t_bwtIndex uiBwtIntervalLength;
		nucSeqIndex uiStartIndexOnQuerry;
		nucSeqIndex uiEndIndexOnQuerry;
		bool bForwHit;
	};

	/* the starting index of the interval (included) */
	nucSeqIndex uiStartIndex;
	/* the end index of the interval (included)*/
	nucSeqIndex uiEndIndex;
	/* list of the perfect matches found through backwards / forward extension */
	std::list<PerfectMatch> lxBwtintervals;
	std::list<PerfectMatch> lxBwtAnchorintervals;
	/* list of the longest perfect matches found through backwards / forward extension */
	std::list<PerfectMatch> lxBwtintervalAnchors;

public:
	SegmentTreeInterval(const nucSeqIndex uiStartIndex, const nucSeqIndex uiEndIndex,
		std::list<PerfectMatch> lxBwtintervals)
		:
		uiStartIndex(uiStartIndex),
		uiEndIndex(uiEndIndex),
		lxBwtintervals(lxBwtintervals)
	{}//constructor
	SegmentTreeInterval(const nucSeqIndex uiStartIndex, const nucSeqIndex uiEndIndex)
		:
		uiStartIndex(uiStartIndex),
		uiEndIndex(uiEndIndex),
		lxBwtintervals()
	{}//constructor
	
	ContainerType getType(){return ContainerType::segment;}//function

	nucSeqIndex getStartIndex() const { return uiStartIndex; }//function
	nucSeqIndex getEndIndex() const { return uiEndIndex; }//function


	/* prints information about this node; thread save */
	void print(std::ostream& xOs) const
	{
		xOs << "(" << std::to_string(this->getStartIndex()) << "," << std::to_string(this->getEndIndex()) << ")";
	}//function
	/* push back an interval of perfect matches
	 * the interval contains uiLengthInBwt individual perfect matches of (uiStartOfIntervalOnQuery, uiEndOfIntervalOnQuery) on the reference sequence
	 * bForwHit is required because we need to extract the starting positions of every single match in the interval using the same fm_index used to calculate the interval
	 * since we use 2 fm_indecies one for forward one for reverse we need to remember which of the two was used
	*/
	void pushBackBwtInterval(t_bwtIndex uiPosInBwt, t_bwtIndex uiLengthInBwt, nucSeqIndex uiStartOfIntervalOnQuery, nucSeqIndex uiEndOfIntervalOnQuery, bool bForwHit, bool bAnchor);

	/* calculates the length of the interval */
	nucSeqIndex length() const
	{
		return uiEndIndex - uiStartIndex + 1;
	}//function

	nucSeqIndex getCenter() const 
	{
		return (uiStartIndex + uiEndIndex + 1) / 2; 
	}//function

	void setInterval(nucSeqIndex uiStart, nucSeqIndex uiEnd) 
	{
		uiStartIndex = uiStart; 
		uiEndIndex = uiEnd; 
	}//function

	unsigned int getValue() const 
	{
		return (unsigned int)length(); 
	}//function

	/* calls fDo for all recorded hits.
	*  Note that pushBackBwtInterval records an interval of hits
	*  ulHit: the position of the hit on the reference sequence
	*  uiQueryBegin: the position of the hit on the query sequence
	*  uiQueryEnd: the end of the hit on the query sequence
	*  bForwHit: true if this hit was produced through forward extension
	*/
	void forEachHitOnTheRefSeq(std::shared_ptr<FM_Index> pxFM_Index, std::shared_ptr<FM_Index> pxRev_FM_Index, 
		unsigned int uiMaxNumHitsPerInterval, bool bSkipLongerIntervals, bool bAnchorOnly,
#if confGENEREATE_ALIGNMENT_QUALITY_OUTPUT
		AlignmentQuality *pxQuality,
#endif
		std::function<void(nucSeqIndex ulIndexOnRefSeq, nucSeqIndex uiQueryBegin, nucSeqIndex uiQueryEnd)> fDo
#ifdef DEBUG_CHECK_INTERVALS
		, NucleotideSequence *pxQuerySeq,
		const BWACompatiblePackedNucleotideSequencesCollection *pxRefSequence
#endif // DEBUG_CHECK_INTERVALS
		)
	{
		//iterate over all the intervals that have been recorded using pushBackBwtInterval()
		for (auto xCurrBwtInterval : bAnchorOnly ? lxBwtAnchorintervals : lxBwtintervals)
		{
			//if the interval contains more than uiMaxNumHitsPerInterval hits it's of no importance and will produce nothing but noise
			t_bwtIndex uiJumpBy = 1;
			if (xCurrBwtInterval.uiBwtIntervalLength > uiMaxNumHitsPerInterval)
			{
				//BOOST_LOG_TRIVIAL(info) << "skipping " << xCurrBwtInterval.uiBwtIntervalLength << " hits for the interval " << *this;
#if confGENEREATE_ALIGNMENT_QUALITY_OUTPUT
				pxQuality->uiAmountOfSkippedBWAIntervals++;
#endif
				if (bSkipLongerIntervals)
					continue;
				uiJumpBy = xCurrBwtInterval.uiBwtIntervalLength / uiMaxNumHitsPerInterval; 
			}//if

			//if the hit was generated using the reversed fm_index we should use the according fm_index in order to 
			//extract the index of the hit on the reference sequence. same for the forward fm_index
			std::shared_ptr<FM_Index> pxUsedFmIndex;
			if (xCurrBwtInterval.bForwHit)
				pxUsedFmIndex = pxRev_FM_Index;
			else
				pxUsedFmIndex = pxFM_Index;

#ifdef DEBUG_CHECK_INTERVALS
			//when debugging check every interval and make sure the matches are correct
			{
				std::unique_lock<std::mutex> xGuard(xMutexDebugCheckIntervals);
				BOOST_LOG_TRIVIAL(info) << "DEBUG_CHECK_INTERVAL on " << *this << " for " << xCurrBwtInterval.uiBwtIntervalLength << " hits bwtIndex=" << xCurrBwtInterval.uiBwtIntervalIndex << (xCurrBwtInterval.bForwHit ? " forw" : " rev");
				SA_IndexInterval ik;
				ik.x[0] = xCurrBwtInterval.uiBwtIntervalIndex;
				ik.x[2] = xCurrBwtInterval.uiBwtIntervalLength;
				ik.setMatchInQuery(xCurrBwtInterval.uiStartIndexOnQuerry, xCurrBwtInterval.uiEndIndexOnQuerry);

				BOOST_LOG_TRIVIAL(info) << ik.startIndexOfMatchInQuery() << " - " << ik.startIndexOfMatchInQuery() + ik.sizeOfMatch() 
					<< " of " << pxQuerySeq->length();

				if (!pxUsedFmIndex->debugCheckInterval(ik, *pxQuerySeq, *pxRefSequence, xCurrBwtInterval.bForwHit))
				{
					BOOST_LOG_TRIVIAL(error) << "DEBUG_CHECK_INTERVAL: Error in: " << *this;
					throw "DEBUG_CHECK_INTERVAL: Error";
				}//if
				else
				{
					BOOST_LOG_TRIVIAL(info) << "DEBUG_CHECK_INTERVAL: aligned " << xCurrBwtInterval.uiBwtIntervalLength 
						<< " hits for the interval " << *this << " correctly";
				}//else
			}//scope of xGuard
#endif
			//iterate over the interval in the BWT
			for (auto ulCurrPos = xCurrBwtInterval.uiBwtIntervalIndex; ulCurrPos < xCurrBwtInterval.uiBwtIntervalIndex + xCurrBwtInterval.uiBwtIntervalLength; ulCurrPos += uiJumpBy)
			{
				//calculate the referenceIndex using pxUsedFmIndex->bwt_sa() and call fDo for every match individually
				auto ulIndexOnRefSeq = pxUsedFmIndex->bwt_sa(ulCurrPos);
				/* if the match was calculated using the fm-index of the reversed sequence:
				 * we acquire the index of the beginning of the match on the reversed sequence by calling bwt_sa()
				 * but we actually want the beginning of the match on the normal sequence, so we need to subtract the END of the match from the reference sequence length
				 */
				if (xCurrBwtInterval.bForwHit)
					ulIndexOnRefSeq = pxUsedFmIndex->getRefSeqLength() - (ulIndexOnRefSeq + xCurrBwtInterval.uiEndIndexOnQuerry - xCurrBwtInterval.uiStartIndexOnQuerry);

				//call the given function
				fDo(ulIndexOnRefSeq, xCurrBwtInterval.uiStartIndexOnQuerry, xCurrBwtInterval.uiEndIndexOnQuerry);
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
		forEachHitOnTheRefSeq(
			pxFM_Index,
			pxRev_FM_Index,
			100000,
			false,
			false,
			[&](nucSeqIndex ulIndexOnRefSeq, nucSeqIndex uiQueryBegin, nucSeqIndex uiQueryEnd)
			{
				vpRet.push_back(pxRefPack->vExtract(
						ulIndexOnRefSeq, 
						ulIndexOnRefSeq + uiQueryEnd - uiQueryBegin - 1
					));
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