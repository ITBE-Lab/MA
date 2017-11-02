/** 
 * @file intervalTree.h
 * @brief Implements a the IntervalTree used for segmentation and various other related classes.
 * @author Markus Schmidt
 */
#ifndef INTERVALTREE_H
#define INTERVALTREE_H

#include "fm_index.h"
#include <thread>
#include "doublyLinkedList.h"
#include "seed.h"

#define confMETA_MEASURE_DURATION ( 1 )


/**
 * @brief A Suffix Array Segment.
 * @details
 * A Suffix Array Segment is made up of two Intervals.
 * @li @c a SA_IndexInterval.
 * @li @c a Interval representing the position of the sequence on the query.
 * @ingroup container
 */
class SaSegment: public Container, public Interval<nucSeqIndex> {
private:
	SA_IndexInterval xSaInterval;
public:
	/**
	* @brief Creates a new SaSegment.
	* @details Creates a new SaSegment on the base of a SA_IndexInterval and the 
	* respective indices on the quey.
	*/
	SaSegment(nucSeqIndex uiStart, nucSeqIndex uiSize, SA_IndexInterval xSaInterval)
			:
		Interval(uiStart, uiSize),
		xSaInterval(xSaInterval)
	{}//constructor

	//overload
	ContainerType getType() const {return ContainerType::segment;}
	/**
	 * @brief The bwt interval within.
	 * @returns the bwt interval within.
	 */
	const SA_IndexInterval& saInterval() const
	{
		return xSaInterval;
	}//function
	
    /**
     * @brief Copys from another SaSegment.
	 * @note override
     */
    inline SaSegment& operator=(const SaSegment& rxOther)
    {
        Interval::operator=(rxOther);
        xSaInterval = rxOther.xSaInterval;
        return *this;
    }// operator

}; // class ( SaSegment )

/**
 * @brief A Interval in the Segment Tree.
 * @ingroup container
 */
class SegmentTreeInterval: public Container, public Interval<nucSeqIndex>
{
private:
	/** 
	 * @brief list of the perfect matches found through backwards / forward extension 
	 */
	std::list<SaSegment> lxSaSegment;

public:
	/**
	 * @brief Creates a new interval with a start and size.
	 */
	SegmentTreeInterval(const nucSeqIndex uiStart, const nucSeqIndex uiSize)
		:
		Interval(uiStart, uiSize),
		lxSaSegment()
	{}//constructor
	
	//overload
	ContainerType getType(){return ContainerType::segment;}//function


	/**
	 * @brief Prints information about this node.
	 * @note Thread save.
	 */
	void print(std::ostream& xOs) const
	{
		xOs << "(" << std::to_string(this->start()) << "," << std::to_string(this->end()) << ")";
	}//function
	/**
	 * @brief Push back an interval of perfect matches.
	 * @details
	 * The interval contains uiLengthInBwt individual perfect matches of 
	 * (uiStartOfIntervalOnQuery, uiEndOfIntervalOnQuery) on the reference sequence.
	 */
	void push_back(SaSegment interval);

	/**
	 * @brief The center of the segment.
	 * @returns the center of the segment.
	 */
	nucSeqIndex getCenter() const 
	{
		return start() + size() / 2; 
	}//function

	/**
	 * @brief Extracts all seeds from the tree.
	 * @details
	 * Calls fDo for all recorded hits.
	 * @Note pushBackBwtInterval records an interval of hits
	 */
	void forEachSeed(
			std::shared_ptr<FM_Index> pxFM_Index,
			unsigned int uiMaxNumHitsPerInterval,
			bool bSkipLongerIntervals,
			std::function<void(Seed s)> fDo
		)
	{
		//iterate over all the intervals that have been recorded using pushBackBwtInterval()
		for (SaSegment xSegment : lxSaSegment)
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

			//iterate over the interval in the BWT
			for (
					auto ulCurrPos = xSegment.saInterval().start(); 
					ulCurrPos < xSegment.saInterval().end(); 
					ulCurrPos += uiJumpBy
				)
			{
				//calculate the referenceIndex using pxUsedFmIndex->bwt_sa() and call fDo for every match individually
				nucSeqIndex ulIndexOnRefSeq = pxFM_Index->bwt_sa(ulCurrPos);

				//TODO: check if this is correct
				//check bridging:
				if(
						ulIndexOnRefSeq*2 < pxFM_Index->getRefSeqLength() && 
						(ulIndexOnRefSeq + xSegment.size() + 1)*2 >= pxFM_Index->getRefSeqLength()
					)
					continue;

				/* if the match was calculated using the fm-index of the reversed sequence:
				 * we acquire the index of the beginning of the match on the reversed sequence by calling bwt_sa()
				 * but we actually want the beginning of the match on the normal sequence, so we need to subtract the END of the match from the reference sequence length
				 */
				//TODO: do i need to throw that away?
				//ulIndexOnRefSeq = pxFM_Index->getRefSeqLength() - (ulIndexOnRefSeq + xSegment.size()) - 1;
				assert(xSegment.start() < xSegment.end());
				//call the given function
				fDo(Seed(xSegment.start(), xSegment.size() + 1, ulIndexOnRefSeq));
			}//for
		}//for
	}//function

	
	/**
	 * @brief returns the number of seeds
	 */
	unsigned int numSeeds()
	{
		return lxSaSegment.size();
	}//function


	/**
	 * @brief Extracts all seeds from the segment.
	 */
	std::shared_ptr<Seeds> getSeeds(
			std::shared_ptr<FM_Index> pxFM_Index, 
			unsigned int min_length
		)
	{
		std::shared_ptr<Seeds> pRet = std::shared_ptr<Seeds>(new Seeds());
		forEachSeed(
				pxFM_Index,
				1000,//parameter has become irrelevant
				true,//parameter has become irrelevant
				[&](Seed xS)
				{
					if(xS.size() >= min_length)
						pRet->push_back(xS);
				}//lambda
			);//forall
		return pRet;
	}//function

	/**
	 * @brief Returns all seeds from the tree.
	 * @details
	 * As opposed to forEachSeed the seeds get collected and returned in a vector.
	 */
	std::vector<std::shared_ptr<NucleotideSequence>> getRefHits(
			std::shared_ptr<FM_Index> pxFM_Index, 
			std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefPack
		)
	{
		std::vector<std::shared_ptr<NucleotideSequence>> vpRet = 
			std::vector<std::shared_ptr<NucleotideSequence>>();
		forEachSeed(
				pxFM_Index,
				10000,
				false,
				[&](Seed xS)
				{
					vpRet.push_back(pxRefPack->vExtract(xS.start_ref(), xS.end_ref()));
				}//lambda
			);//forall
		return vpRet;
	}//function
};//class

/**
 * @brief The segment tree.
 * @details
 * The segment "tree" is actually a doubly linked list.
 * The tree only exists logically,
 * meaning that the segments within the list represent the first layer of the tree initally.
 * Then after each iteration, the segments within the list represent the next layer down of the 
 * tree.
 * @ingroup container
 */
class SegmentTree : public DoublyLinkedList<SegmentTreeInterval>, public Container{

public:
	/**
	* @brief Creates a new tree containing one initial segment as root.
	* @details
	* Sets up the interval tree with two leaves and one initial interval comprising the whole query
	* note that the tree is internally represented as a DoublyLinkedList since only the leaves are of relevance
	*/
	SegmentTree(const nucSeqIndex uiQueryLength)
		:
		DoublyLinkedList()
	{
		//the intervals are inclusive in the tree...
		std::shared_ptr<SegmentTreeInterval> pxRoot(new SegmentTreeInterval(0, uiQueryLength-1));
		push_back(pxRoot);
	}//constructor
	SegmentTree()
	{}//constructor
	

    static std::shared_ptr<SegmentTree> cast(std::shared_ptr<Container> p)
    {
        return std::dynamic_pointer_cast<SegmentTree>(p);
    }
	
	//overload
	ContainerType getType(){return ContainerType::segmentList;}//function

	/**
	 * @brief Prints basic information about the segment tree.
	 * @details
	 * Not thread save.
	 */
	void print(std::ostream &xOut) const
	{
		forEach([&xOut](std::shared_ptr<SegmentTreeInterval> pxNode){ pxNode->print(xOut); });
	}//function

	
	/**
	 * @brief Extracts all seeds from the segment.
	 */
	std::shared_ptr<Seeds> getSeeds(
			std::shared_ptr<FM_Index> pxFM_Index, 
			unsigned int min_length
		)
	{
		
		std::shared_ptr<Seeds> pRet = std::shared_ptr<Seeds>(new Seeds());
		forEach(
			[&]
			(std::shared_ptr<SegmentTreeInterval> pxNode)
			{ 
				pRet->append(pxNode->getSeeds(pxFM_Index, min_length));
			}//lambda
		);
		return pRet;
	}//function

	
	/**
	 * @brief returns the number of seeds
	 */
	unsigned int numSeeds()
	{
		
		unsigned int uiTotal = 0;
		forEach(
			[&]
			(std::shared_ptr<SegmentTreeInterval> pxNode)
			{ 
				uiTotal += pxNode->numSeeds();
			}//lambda
		);
		return uiTotal;
	}//function
};

/**
 * @brief Simple printer function for the SegmentTree.
 */
std::ostream& operator<<(std::ostream& xOs, const SegmentTree& rxTree);
/**
 * @brief Simple printer function for a SegmentTreeInterval.
 */
std::ostream& operator<<(std::ostream& xOs, const SegmentTreeInterval &rxNode);

/**
 * @brief Exposes the SegmentTree to boost python.
 * @ingroup export
 */
void exportIntervalTree();


#endif