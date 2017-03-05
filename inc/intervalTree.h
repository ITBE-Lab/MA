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

//#define DEBUG_ALIGNER 
//#define DEBUG_CHECK_INTERVALS 
//changing this to 1 might break a lot of stuff (old code) [i don't think it's worth it to go through and repair it though]
#define USE_BUCKET_CLUTERING ( 0 )

#define confMETA_MEASURE_DURATION ( 1 )
#define confGENEREATE_ALIGNMENT_QUALITY_OUTPUT ( 1 )

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

template <typename Content>
class DoublyLinkedList;

/* The elements of the doubly linked list.
*/
template <typename Content>
class DoublyLinkedListElement
{
private:
	/* pointer to the the last element of the linked list */
	std::weak_ptr<DoublyLinkedListElement<Content>> pxLastNode;
	/* the content of the list element */
	std::shared_ptr<Content> pxContent;
	/* pointer to the next element of the linked list */
	std::shared_ptr<DoublyLinkedListElement<Content>> pxNextNode;

	//Deadlock prevention technique: always lock the elements IN-order of the list
	/* when multi threading the nodes need to be locked before inserting removing elements of the list */
	std::mutex xMutex;

public:
	DoublyLinkedListElement(std::weak_ptr<DoublyLinkedListElement<Content>> pxLastNode,
		std::shared_ptr<Content> pxContent, std::shared_ptr<DoublyLinkedListElement<Content>> pxNextNode)
		:
		pxLastNode(pxLastNode),
		pxContent(pxContent),
		pxNextNode(pxNextNode)
	{}//constructor
	DoublyLinkedListElement(){}//constructor

	//not thread save!
	void setNextNode(std::shared_ptr<DoublyLinkedListElement<Content>> pxNextNode)
	{
		this->pxNextNode = pxNextNode;
	}//function
	//not thread save!
	void setLastNode(std::weak_ptr<DoublyLinkedListElement<Content>> pxLastNode)
	{
		this->pxLastNode = pxLastNode;
	}//function
	const std::shared_ptr<DoublyLinkedListElement<Content>> getNextNode() const { return pxNextNode; }//function
	const std::shared_ptr<DoublyLinkedListElement<Content>> getLastNode() const { return pxLastNode.lock(); }//function
	std::shared_ptr<DoublyLinkedListElement<Content>> getNextNode() { return pxNextNode; }//function
	std::shared_ptr<DoublyLinkedListElement<Content>> getLastNode() { return pxLastNode.lock(); }//function
	/* the user of the linked list will never have access to an object of DoublyLinkedListElement so the mutex will be used from DoublyLinkedList only */
	std::mutex *getMutex() { return &xMutex; }//function
	/* the two ends of the doubly linked list are object from the class DoublyLinkedListEnd; 
	 * this method is used to distinguish between DoublyLinkedListEnd and DoublyLinkedListElement 
	 */
	virtual bool isListElement() const { return true; }//function
	virtual std::shared_ptr<Content> getContent() const { return pxContent; }//function
};//class

template <typename Content>
class DoublyLinkedListEnd : public DoublyLinkedListElement<Content>
{

public:
	DoublyLinkedListEnd()
		:
		DoublyLinkedListElement<Content>()
	{}

	bool isListElement() const { return false; };//function
	std::shared_ptr<Content> getContent() const { return nullptr; }//function
};//class

template <typename Content>
class DoublyLinkedList
{

private:
	/* points to the right leaf of the list; used to determine the end of the list. */
	std::shared_ptr<DoublyLinkedListEnd<Content>> pxLastLeaf = nullptr;
	/* points to the left leaf of the list; used to determine the beginning of the list.*/
	std::shared_ptr<DoublyLinkedListEnd<Content>> pxFrirstLeaf = nullptr;

	/* returns the first element (containing data) of the list;*/
	const std::shared_ptr<DoublyLinkedListElement<Content>> getRoot() const { return pxFrirstLeaf->getNextNode(); }
	/* returns the first element (containing data) of the list;*/
	std::shared_ptr<DoublyLinkedListElement<Content>> getRoot() { return pxFrirstLeaf->getNextNode(); }
public:

	/* objects of Iterator are used to give any users of the list access to the individual list elements
	*/
	class Iterator{
		/* 
		 * allow doubly linked list to access all private members
		 * necessary because we need to access the actual DoublyLinkedListElement in DoublyLinkedList
		 * but all other classes should not have access to anything but the Content
		*/
		friend class DoublyLinkedList;
	private:
		std::shared_ptr<DoublyLinkedListElement<Content>> pxPos;

		std::shared_ptr<DoublyLinkedListElement<Content>> getElement(){ return pxPos; }//function
	public:
		Iterator(std::shared_ptr<DoublyLinkedListElement<Content>> pxPos)
			:
			pxPos(pxPos)
		{}//constructor
		Iterator(const Iterator &rxOther)
			:
			pxPos(rxOther.pxPos)
		{}//constructor

		void operator++(){ pxPos = pxPos->getNextNode(); }//function
		void operator--(){ pxPos = pxPos->getLastNode(); }//function
		Iterator getCopy() const { return Iterator(*this); }//function
		void operator=(const Iterator &rxOther){ pxPos = rxOther.pxPos; }//function
		std::shared_ptr<Content> operator*() const { return pxPos->getContent(); }//function
		Content* operator->() const { return pxPos->getContent().get(); }//function
		bool isListElement() const { return pxPos->isListElement(); }//function
	};
	/*
	*	sets up an empty doubly linked list
	*/
	DoublyLinkedList()
		:
		pxFrirstLeaf(new DoublyLinkedListEnd<Content>()),
		pxLastLeaf(new DoublyLinkedListEnd<Content>())
	{
		pxFrirstLeaf->setNextNode(pxLastLeaf);
		pxLastLeaf->setLastNode(pxFrirstLeaf);
	}//constructor

	Iterator begin() const { return Iterator(getRoot()); }//function
	Iterator end() const { return Iterator(pxLastLeaf->getLastNode()); }//function

	void forEach(std::function<void(std::shared_ptr<Content>)> fDo) const
	{
		std::shared_ptr<DoublyLinkedListElement<Content>> pxCurr = getRoot();
		while (pxCurr->isListElement())
		{
			fDo(pxCurr->getContent());
			pxCurr = pxCurr->getNextNode();
		}//while
	}//function

	void removeNode(Iterator pxItt)
	{
		removeNode(pxItt.getElement());
	}//function

	void removeNode(std::shared_ptr<DoublyLinkedListElement<Content>> pxNode)
	{
		//lock the last, this and the next node
		//Deadlock prevention technique: always lock the elements IN-order of the list
		std::lock_guard<std::mutex> xOtherGuard(*pxNode->getLastNode()->getMutex());
		std::lock_guard<std::mutex> xMyGuard(*pxNode->getMutex());
		std::lock_guard<std::mutex> xOtherGuard2(*pxNode->getNextNode()->getMutex());
		//check if this node is the root node; if so move the root node to the next node
		pxNode->getLastNode()->setNextNode(pxNode->getNextNode());
		pxNode->getNextNode()->setLastNode(pxNode->getLastNode());
	}//function

	Iterator insertBefore(std::shared_ptr<Content> pxContent, std::shared_ptr<DoublyLinkedListElement<Content>> pxNode)
	{
		//lock the last and this node
		//Deadlock prevention technique: always lock the elements IN-order of the list
		std::lock_guard<std::mutex> xOtherGuard(*pxNode->getLastNode()->getMutex());
		std::lock_guard<std::mutex> xMyGuard(*pxNode->getMutex());

		std::shared_ptr<DoublyLinkedListElement<Content>> pxNewElement
			= std::shared_ptr<DoublyLinkedListElement<Content>>(new DoublyLinkedListElement<Content>(pxNode->getLastNode(), pxContent, pxNode));

		pxNode->getLastNode()->setNextNode(pxNewElement);
		pxNode->setLastNode(pxNewElement);

		return Iterator(pxNewElement);
	}//function

	Iterator insertBefore(std::shared_ptr<Content> pxContent, Iterator pxItt)
	{
		return insertBefore(pxContent, pxItt.getElement());
	}//function

	Iterator insertAfter(std::shared_ptr<Content> pxContent, std::shared_ptr<DoublyLinkedListElement<Content>> pxNode)
	{
		//lock the next and this node
		//Deadlock prevention technique: always lock the elements IN-order of the list
		std::lock_guard<std::mutex> xMyGuard(*pxNode->getMutex());
		std::lock_guard<std::mutex> xOtherGuard2(*pxNode->getNextNode()->getMutex());

		std::shared_ptr<DoublyLinkedListElement<Content>> pxNewElement
			= std::shared_ptr<DoublyLinkedListElement<Content>>(new DoublyLinkedListElement<Content>(pxNode, pxContent, pxNode->getNextNode()));

		pxNode->getNextNode()->setLastNode(pxNewElement);
		pxNode->setNextNode(pxNewElement);
			
		return Iterator(pxNewElement);
	}//function

	Iterator insertAfter(std::shared_ptr<Content> pxContent, Iterator pxItt)
	{
		return insertAfter(pxContent, pxItt.getElement());
	}//function

	Iterator push_back(std::shared_ptr<Content> pxContent)
	{
		return insertBefore(pxContent, pxLastLeaf);
	}//function
	Iterator push_front(std::shared_ptr<Content> pxContent)
	{
		return insertAfter(pxContent, pxFrirstLeaf);
	}//function

	/*the number of elements in the list; non thread save*/
	unsigned int length() const
	{
		unsigned int uiRet = 0;

		forEach([&uiRet](std::shared_ptr<Content>){ uiRet++; });//forEach

		return uiRet;
	}//function

};//class

class SegmentTreeInterval{
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
	void pushBackBwtInterval(t_bwtIndex uiPosInBwt, t_bwtIndex uiLengthInBwt, nucSeqIndex uiStartOfIntervalOnQuery, nucSeqIndex uiEndOfIntervalOnQuery, bool bForwHit, bool bAnchor)
	{
		PerfectMatch xMatch =
		{
			uiPosInBwt, // uiBwtIntervalIndex
			uiLengthInBwt, // uiBwtIntervalLength
			uiStartOfIntervalOnQuery, // uiStartIndexOnQuerry
			uiEndOfIntervalOnQuery, // uiEndIndexOnQuerry
			bForwHit, // bForwHit
		};
		lxBwtintervals.push_back(xMatch);
		if (bAnchor)
			lxBwtAnchorintervals.push_back(xMatch);
	}//function

	/* calculates the length of the interval */
	nucSeqIndex length() const
	{
		return uiEndIndex - uiStartIndex + 1;
	}//function
	nucSeqIndex getCenter() const { return (uiStartIndex + uiEndIndex + 1) / 2; }//function
	void setInterval(nucSeqIndex uiStart, nucSeqIndex uiEnd) { uiStartIndex = uiStart; uiEndIndex = uiEnd; }//function
	unsigned int getValue() const { return (unsigned int)length(); }//function
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
};

class SegmentTree : public DoublyLinkedList<SegmentTreeInterval>{

public:
	/*
	*	sets up the interval tree with two leaves and one initial interval comprising the whole query
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

	/*not thread save; prints basic information about the segment tree*/
	void print(std::ostream &xOut) const
	{
		forEach([&xOut](std::shared_ptr<SegmentTreeInterval> pxNode){ pxNode->print(xOut); });
	}//function


	SegmentTree getTheNLongestIntervals(unsigned int uiN)
	{
		SegmentTree xRet;

		forEach(
			[&xRet, &uiN](std::shared_ptr<SegmentTreeInterval> pxNode)
			{
				auto pxIterator = xRet.begin();
				while (pxIterator.isListElement() && pxIterator->length() > pxNode->length())
					++pxIterator;
				xRet.insertBefore(pxNode, pxIterator);
				if (xRet.length() > uiN)
					xRet.removeNode(xRet.end());
			}//lambda
		);//forEach

		return xRet;
	}//function
};

std::ostream& operator<<(std::ostream& xOs, const SegmentTree& rxTree);
std::ostream& operator<<(std::ostream& xOs, const SegmentTreeInterval &rxNode);

#if 0

/*
*	A node of the interval tree;
*	the tree is internally represented as a linked list,
*		since only the leaves of the tree are of any relevance
*		and the tree structure contains no valuable information
*	One node represents one interval.
*	class is thread save
*/
class SegmentTreeNode {

protected:
	/* the starting index of the interval (included) */
	unsigned int uiStartIndex;
	/* the end index of the interval (included)*/
	unsigned int uiEndIndex;
	/* list of the hits from backwards extension in this interval.
	*	tuple: (index in bwt, length in bwt, start pos of the interval on the query) 
	*	[the end pos of an interval is never changed]
	*/
	std::list<std::tuple<t_bwtIndex, t_bwtIndex, unsigned int>> ltuibwtInterval;


#ifdef DEBUG_INTERVALTREE
	bool bBwtValsSet = false;
#endif //DEBUG_INTERVALTREE
	unsigned int uiBackwardsExtensionReachedUntil = 0;
	bool bHasStoredBackwardsExtension = false;
	/* pointers to the next and the last element of the linked list */
	std::shared_ptr<SegmentTreeNode> pxNextNode;
	std::weak_ptr<SegmentTreeNode> pxLastNode;
	/* when multi threading the nodes need to be locked before splitting */
	std::mutex xMutex;

#ifdef DEBUG_INTERVALTREE
	/* when debugging we give each element an unique id in order to track it */
	unsigned int uiId = iNextId++;
#endif // DEBUG_INTERVALTREE

#ifdef DEBUG_ALIGNER
	//just used for debugging
	bool bHitsSaved = false;
#endif // DEBUG_ALIGNER

public:
	SegmentTreeNode(const unsigned int uiStartIndex, const unsigned int uiEndIndex, unsigned int uiLastError, 
		std::shared_ptr<SegmentTreeNode> pxNextNode, std::weak_ptr<SegmentTreeNode> pxLastNode);

	SegmentTreeNode(const unsigned int uiStartIndex, const unsigned int uiEndIndex,
		std::shared_ptr<SegmentTreeNode> pxNextNode, std::weak_ptr<SegmentTreeNode> pxLastNode);

	SegmentTreeNode(const unsigned int uiStartIndex, const unsigned int uiEndIndex,
		std::shared_ptr<SegmentTreeNode> pxNextNode, std::weak_ptr<SegmentTreeNode> pxLastNode,
		std::list<std::tuple<t_bwtIndex, t_bwtIndex, unsigned int>>);


	SegmentTreeNode(const unsigned int uiStartIndex, const unsigned int uiEndIndex);

	SegmentTreeNode(); 
	~SegmentTreeNode();

	/* calls print on this and all following nodes; not thread save*/
	virtual void printTree(std::ostream& xOs) const;
	/* prints information about this node; thread save */
	virtual void print(std::ostream& xOs) const;
	/* calculates the length of the interval; thread save */
	virtual unsigned int length() const;
	/* splits the interval exactly in half; in order to do this one new node is added to the linked list and returned for confinence; is thread save*/
	virtual std::shared_ptr<SegmentTreeNode> split(unsigned int uiSplitIndex, std::shared_ptr<SegmentTreeNode> pxThis);

	virtual unsigned int getValue() const { return length(); }//function
	/* removes this node from the list; threadsave*/
	virtual void removeNode(SegmentTree &rxTree, std::shared_ptr<SegmentTreeNode> pxThis);
	//setter not thread save!
	void setNextNode(std::shared_ptr<SegmentTreeNode> pxNextNode);
	void setLastNode(std::weak_ptr<SegmentTreeNode> pxLastNode);
	void setInterval(unsigned int uiStart, unsigned int uiEnd) { uiStartIndex = uiStart; uiEndIndex = uiEnd; }//function
	void forEachHitOnTheRefSeq(std::shared_ptr<FM_Index> pxFM_Index, std::function<void(unsigned long ulHit, unsigned int uiQueryIndex)> fDo
#ifdef DEBUG_CHECK_INTERVALS
		, const NucleotideSequence *pxQuerySeq, const BWACompatiblePackedNucleotideSequencesCollection *pxRefSequence
#endif // DEBUG_CHECK_INTERVALS
		);
	//getter not thread save!
	const unsigned int getStartIndex() const { return uiStartIndex; }//function
	const unsigned int getEndIndex() const { return uiEndIndex; }//function
	const std::shared_ptr<SegmentTreeNode> getNextNode() const { return pxNextNode; }//function
	const std::shared_ptr<SegmentTreeNode> getLastNode() const { return pxLastNode.lock(); }//function
	std::shared_ptr<SegmentTreeNode> getNextNode() { return pxNextNode; }//function
	std::shared_ptr<SegmentTreeNode> getLastNode() { return pxLastNode.lock(); }//function
	void pushBackBwtInterval(t_bwtIndex uiPosInBwt, t_bwtIndex uiLengthInBwt, unsigned int uiStartOfIntervalOnQuery);
	unsigned int getCenter() const { return (uiStartIndex + uiEndIndex + 1) / 2; }//function

	void storeBackwardsExtension(unsigned int uiError) { bHasStoredBackwardsExtension = true; this->uiBackwardsExtensionReachedUntil = uiError; }//function
	bool hasStoredBackwardsExtension() const { return bHasStoredBackwardsExtension; }//function
	unsigned int backwardsExtensionReachedUntil() { return uiBackwardsExtensionReachedUntil; }//function
	virtual bool isNode() { return true; }//function

#ifdef DEBUG_ALIGNER
	bool getHitsSaved() const { return bHitsSaved; }//function
	void setHitsSaved() { bHitsSaved = true; }//function
#endif //DEBUG_ALIGNER
};//class
/*
*	one "leaf" at each end of the list;
*	the leafs do not contain any data
*	used for the termination of recursive calls
*	the next pointer of the last leaf points to itself 
*	the last pointer of the first leaf points to itself
*/
class SegmentTreeLeaf : public SegmentTreeNode {

public:
	SegmentTreeLeaf()
		:
		SegmentTreeNode()
	{}//constructor

	SegmentTreeLeaf(std::shared_ptr<SegmentTreeNode> pxNextNode, std::weak_ptr<SegmentTreeNode> pxLastNode)
			:
		SegmentTreeNode(0, 0, pxNextNode, pxLastNode)
	{}//constructor

	void printTree(std::ostream& xOs) const {}//function
	void print(std::ostream& xOs) const {}//function
	unsigned int getValue() const { return 0; }//function

	unsigned int length() const { return 0; }//function

	std::shared_ptr<SegmentTreeNode> split(unsigned int uiSplitIndex, std::shared_ptr<SegmentTreeNode> pxThis) { return nullptr; }//function
	void removeNode() {}//function
	bool isNode() { return false; }//function
};
#endif

#endif