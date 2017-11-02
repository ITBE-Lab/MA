#if 0
#include "reSeed.h"


#include <vector>
#include <memory>
#include <atomic>
#include <chrono>
 
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
}//function


std::vector<ContainerType> ReSeed::getInputType()
{
	return std::vector<ContainerType>{
			//the forward fm_index
			ContainerType::fM_index,
			//the forward fm_index
			ContainerType::segmentList,
			//the query sequence
			ContainerType::nucSeq,
		};
}
ContainerType ReSeed::getOutputType()
{
	return ContainerType::segmentList;
}


std::shared_ptr<Container> ReSeed::execute(
		std::vector<std::shared_ptr<Container>> vpInput
	)
{
	std::shared_ptr<FM_Index> pFM_index = std::static_pointer_cast<FM_Index>(vpInput[0]);
	std::shared_ptr<SegmentTree> pSegments = std::static_pointer_cast<SegmentTree>(vpInput[1]);
	std::shared_ptr<NucleotideSequence> pQuerySeq = 
		std::static_pointer_cast<NucleotideSequence>(vpInput[2]);

		
	std::shared_ptr<SegmentTree> pSegmentTree(new SegmentTree(pQuerySeq->length()));



	return pSegmentTree;
}//function

void exportReSeed()
{
	//export the ReSeed class
	boost::python::class_<ReSeed, boost::python::bases<CppModule>>("ReSeed")
		.read_write("min_split_len", &ReSeed::minSplitLen)
	;
}//function
#endif