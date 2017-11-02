#include "longestNonEnclosedSegments.h"

#define INCLUDE_SELF_TESTS (  1 )

//#include "BWTmem.h"
#include <vector>
//#include "analysis.h"
#include <memory>
#include <atomic>
#include <chrono>
//#include "assembly.h"
 
#define complement(x) (uint8_t)NucleotideSequence::nucleotideComplement(x)

/**
* Delivers 4 intervals for a single input interval.
* Here we use only two fields in the BWT_Interval.
*/
SA_IndexInterval LongestNonEnclosedSegments::extend_backward( 
		// current interval
		const SA_IndexInterval &ik,
		// the character to extend with
		const uint8_t c,
		std::shared_ptr<FM_Index> pFM_index
	)
{
	bwt64bitCounter cntk[4]; // Number of A, C, G, T in BWT until start of interval ik
	bwt64bitCounter cntl[4]; // Number of A, C, G, T in BWT until end of interval ik

	assert(ik.start() < ik.end());
	assert(ik.start() > 0);

	//here the intervals seem to be (a,b] while mine are [a,b)
	pFM_index->bwt_2occ4(
		// start of SA index interval
		ik.start() - 1,
		// end of SA index interval
		ik.end() - 1,
		cntk,						// output: Number of A, C, G, T until start of interval
		cntl						// output: Number of A, C, G, T until end of interval
	);

	for(unsigned int i = 0; i < 4; i++)
		assert(cntk[i] <= cntl[i]);

	bwt64bitCounter cnts[4]; // Number of A, C, G, T in BWT interval ik
	//the cnts calculated here might be off by one
	for(unsigned int i = 0; i < 4; i++)
		cnts[i] = cntl[i] - cntk[i];

	DEBUG_2(
		std::cout << cnts[0] << " + " << cnts[1] << " + " << cnts[2] << " + " << cnts[3] << " = " 
				  << (t_bwtIndex)(cnts[0] + cnts[1] + cnts[2] + cnts[3]) << " ?= " 
				  << ik.size() << "(-1)" << std::endl;
	)

	bwt64bitCounter cntk_2[4];
	cntk_2[0] = ik.startRevComp();
	/*
	 * PROBLEM:
	 * 
	 * the representation of the $ in the count part of the FM_index is indirect
	 * 		done by storing the position of the $
	 * if have two bwt indices k and l
	 * the counts do not return the $ obviously...
	 * 
	 * The result may be off by one since sometimes we have a $ before the current pos 
	 * sometimes we do not...
	 *
	 * lets adjust the sizes of the smaller intervals accordingly
	 */
	if(
			ik.start() < pFM_index->primary && 
			ik.end() >= pFM_index->primary
		)
	{
		cntk_2[0]++;
		DEBUG_2(
			std::cout << "adjusted cntk_2[0] because of primary" << std::endl;
		)
		assert( (t_bwtIndex)(cnts[0] + cnts[1] + cnts[2] + cnts[3]) == ik.size() - 1 );
	}//if
	else{
		if( (t_bwtIndex)(cnts[0] + cnts[1] + cnts[2] + cnts[3]) != ik.size() )
		{
			std::cout << ik.start() << " " << ik.end() << " " << pFM_index->primary << std::endl;
			std::cout << cnts[0] << " + " << cnts[1] << " + " << cnts[2] << " + " <<
				cnts[3] << " = " << (t_bwtIndex)(cnts[0] + cnts[1] + cnts[2] + cnts[3]) << " ?= "
				<< ik.size() << "(-1)" << std::endl;
		}//if
		assert( (t_bwtIndex)(cnts[0] + cnts[1] + cnts[2] + cnts[3]) == ik.size() );
	}//else
	//for all nucleotides
	for(unsigned int i = 1; i < 4; i++)
		cntk_2[i] = cntk_2[i-1] + cnts[complement(i-1)];



	//BWAs SA intervals seem to be (a,b] while mine are [a,b)
	//pFM_index->L2[c] start of nuc c in BWT
	//cntk[c] offset of new interval
	//cntl[c] end of new interval
	return SA_IndexInterval(pFM_index->L2[c] + cntk[c] + 1, cntk_2[complement(c)], cnts[c]);
} // method


SaSegment LongestNonEnclosedSegments::extend(
		std::shared_ptr<SegmentTreeInterval> pxNode,
		std::shared_ptr<FM_Index> pFM_index,
		std::shared_ptr<NucleotideSequence> pQuerySeq
	)
{
	nucSeqIndex center = pxNode->start() + pxNode->size()/2;

	// query sequence itself
	const uint8_t *q = pQuerySeq->pGetSequenceRef(); 
	
	/* Initialize ik on the foundation of the single base q[x].
	 * In order to understand this initialization you should have a look 
	 *to the corresponding PowerPoint slide.
	 */
	// start I(q[x]) in T (start in BWT used for backward search) + 1, 
	// because very first string in SA-array starts with $
	// size in T and T' is equal due to symmetry
	SA_IndexInterval ik(
						pFM_index->L2[complement(q[center])] + 1, 
						pFM_index->L2[(int)q[center]] + 1, 
						pFM_index->L2[(int)q[center] + 1] - pFM_index->L2[(int)q[center]]
					);

	std::list<SaSegment> curr = std::list<SaSegment>();
	for(nucSeqIndex i = center+1; i < pQuerySeq->length(); i++)
	{
		DEBUG_2(
			std::cout << i-1 << " -> " << ik.start() << " " << ik.end() << std::endl;
			std::cout << i-1 << " ~> " << ik.revComp().start() << " " << ik.revComp().end() << std::endl;
		)
		assert(ik.size() > 0);
		SA_IndexInterval ok = extend_backward(ik, complement(q[i]), pFM_index);

		if(ok.size() != ik.size()) 
			curr.push_front(SaSegment(center, i-center-1, ik.revComp()));
		if(i == pQuerySeq->length()-1 && ok.size() != 0)
			curr.push_front(SaSegment(center, i-center, ok.revComp()));

		DEBUG_2(
			std::cout << i << " -> " << ok.start() << " " << ok.end() << std::endl;
			std::cout << i << " ~> " << ok.revComp().start() << " " << ok.revComp().end() << std::endl;
		)
		/*
		* In fact, if ok.getSize is zero, then there are no matches any more.
		*/
		if (ok.size() == 0)
			break; // the SA-index interval size is too small to be extended further
		ik = ok;
	}//for
	DEBUG_2(
		std::cout << "swap" << std::endl;
	)
	std::list<SaSegment> prev = std::list<SaSegment>();
	SaSegment longest(0,0,SA_IndexInterval(0,0,0));
	std::list<SaSegment> *pPrev, *pCurr, *pTemp;
	pPrev = &curr;
	pCurr = &prev;
	if(center != 0)
	{
		for(nucSeqIndex i = center-1; i >= 0; i--)
		{
			assert(pCurr->empty());

			bool bHaveOne = false;

			for(SaSegment ik : *pPrev)
			{
				DEBUG_2(
					std::cout << i+1 << " -> " << ik.saInterval().start() << " " << ik.saInterval().end() << std::endl;
					std::cout << i+1 << " ~> " << ik.saInterval().revComp().start() << " " << ik.saInterval().revComp().end() << std::endl;
				)
				SA_IndexInterval ok = extend_backward(ik.saInterval(), q[i], pFM_index);
				DEBUG_2(
					std::cout << i << " -> " << ok.start() << " " << ok.end() << std::endl;
					std::cout << i << " ~> " << ok.revComp().start() << " " << ok.revComp().end() << std::endl;
				)
				DEBUG(
					std::cout << ik.start() << ", " << ik.end() << ": " << ik.saInterval().size() << " -> " << ok.size() << std::endl;
				)
				if(ok.size() == 0 && !bHaveOne)
				{
					pxNode->push_back(ik);
					assert(ik.end() <= pQuerySeq->length());
					if(ik.size() > longest.size())
						longest = ik;
					bHaveOne = true;
				}//if
				else if(ok.size() != 0)
				{
					SaSegment xSeg = SaSegment(i, ik.size()+1, ok);
					pCurr->push_back(xSeg);
					assert(xSeg.end() <= pQuerySeq->length());
				}//if
			}//for
			pTemp = pPrev;
			pPrev = pCurr;
			pCurr = pTemp;
			pCurr->clear();

			//if there are no more intervals to extend
			if(pPrev->empty())
				break;

			//cause nuxSeqIndex is unsigned
			if(i == 0)
				break;
		}//for
	}//if
	
	if(!pPrev->empty())
	{
		assert(pPrev->front().size() >= pPrev->back().size());

		pxNode->push_back(pPrev->front());
		assert(pPrev->front().end() <= pQuerySeq->length());
		
		DEBUG_2(
			std::cout << pPrev->front().start() << ":" << pPrev->front().end() << std::endl;
		)

		if(pPrev->front().size() > longest.size())
			longest = pPrev->front();
	}//if

	return longest;

	/*if (bBackwards)
		return SaSegment(i + 1, uiStartIndex - (i + 1), ik, false);
	else
		return SaSegment(uiStartIndex, i - uiStartIndex - 1, ik, true);*/
}//function

/* this function implements the segmentation of the query
 *
 * the process is synchronized using a thread pool
 *
 * to avoid unnecessary locking and unlocking we make sure that only one thread can process one interval at any given time.
 * ( there are always several intervals ready and waiting to be processed. )
 * this way it is only necessary to lock once we touch the structure of the list which stores the individual intervals.
 *
 * segmentation technique:
 *		start in the middle of the interval try to extend in both directions
 *		split the interval in 3 parts: prev:perfectMatch:post
 *		queue the prev and post intervals into the thread pool
 *		save the perfect match for later clustering
*/
void LongestNonEnclosedSegments::procesInterval(
			size_t uiThreadId,
			DoublyLinkedList<SegmentTreeInterval>::Iterator pxNode,
			std::shared_ptr<SegmentTree> pSegmentTree,
			std::shared_ptr<FM_Index> pFM_index,
			std::shared_ptr<NucleotideSequence> pQuerySeq,
			ThreadPoolAllowingRecursiveEnqueues* pxPool
		)
{
	DEBUG(
		std::cout << "interval (" << pxNode->start() << "," << pxNode->end() << ")" << std::endl;
	)

	//performs extension and records any perfect matches
	SaSegment xLongest = extend(*pxNode, pFM_index, pQuerySeq);

	nucSeqIndex uiFrom, uiTo;
	uiFrom = xLongest.start();
	uiTo = xLongest.end();
	DEBUG(
		std::cout << "splitting interval (" << pxNode->start() << "," << pxNode->end() << ") at (" << uiFrom << "," << uiTo << ")" << std::endl;
	)

	if (uiFrom != 0 && pxNode->start() + 1 < uiFrom)
	{
		//create a new list element and insert it before the current node
		auto pxPrevNode = pSegmentTree->insertBefore(std::shared_ptr<SegmentTreeInterval>(
			new SegmentTreeInterval(pxNode->start(), uiFrom - pxNode->start() - 1)), pxNode);
		//enqueue procesInterval() for the new interval
		pxPool->enqueue( 
			LongestNonEnclosedSegments::procesInterval,
			pxPrevNode, pSegmentTree, pFM_index, pQuerySeq, pxPool
		);//enqueue
	}//if
	if (pxNode->end() > uiTo + 1)
	{
		//create a new list element and insert it after the current node
		auto pxNextNode = pSegmentTree->insertAfter(std::shared_ptr<SegmentTreeInterval>(
			new SegmentTreeInterval(uiTo + 1, pxNode->end() - uiTo - 1)), pxNode);
		//enqueue procesInterval() for the new interval
		pxPool->enqueue( 
			LongestNonEnclosedSegments::procesInterval,
			pxNextNode, pSegmentTree, pFM_index, pQuerySeq, pxPool
		);//enqueue
	}//if


	pxNode->start(uiFrom);
	pxNode->end(uiTo);
}//function


std::vector<ContainerType> LongestNonEnclosedSegments::getInputType()
{
	return std::vector<ContainerType>{
			//the forward fm_index
			ContainerType::fM_index,
			//the query sequence
			ContainerType::nucSeq,
		};
}
ContainerType LongestNonEnclosedSegments::getOutputType()
{
	return ContainerType::segmentList;
}


std::shared_ptr<Container> LongestNonEnclosedSegments::execute(
		std::vector<std::shared_ptr<Container>> vpInput
	)
{
	std::shared_ptr<FM_Index> pFM_index = std::static_pointer_cast<FM_Index>(vpInput[0]);
	std::shared_ptr<NucleotideSequence> pQuerySeq = 
		std::static_pointer_cast<NucleotideSequence>(vpInput[1]);

		
	std::shared_ptr<SegmentTree> pSegmentTree(new SegmentTree(pQuerySeq->length()));


	assert(*pSegmentTree->begin() != nullptr);

	{//scope for xPool
		DEBUG(
			ThreadPoolAllowingRecursiveEnqueues xPool( 1 );
		)
		#if DEBUG_LEVEL <= 0
			ThreadPoolAllowingRecursiveEnqueues xPool( NUM_THREADS_ALIGNER );
		#endif

		//enqueue the root interval for processing
		xPool.enqueue( 
			LongestNonEnclosedSegments::procesInterval,
			pSegmentTree->begin(), pSegmentTree, pFM_index, pQuerySeq, &xPool
		);//enqueue

	}//end of scope xPool

	return pSegmentTree;
}//function

void exportLongestNonEnclosedSegments()
{
	//boost::python::def("analyse_crisper", analyseCRISPER);



	//export the LongestNonEnclosedSegments class
	boost::python::class_<LongestNonEnclosedSegments, boost::python::bases<CppModule>>(
			"LongestNonEnclosedSegments",
			"bBreakOnAmbiguousBase: weather the extension of "
			"intervals shall be stopped at N's\n"
		)
		;
	
}//function