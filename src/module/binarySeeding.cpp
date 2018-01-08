#include "module/binarySeeding.h"

using namespace libMABS;

#define INCLUDE_SELF_TESTS (  1 )

//#include "BWTmem.h"
#include <vector>
//#include "analysis.h"
#include <memory>
#include <atomic>
#include <chrono>
//#include "assembly.h"
 
#define complement(x) (uint8_t)NucSeq::nucleotideComplement(x)

SAInterval BinarySeeding::extend_backward( 
		// current interval
		const SAInterval &ik,
		// the character to extend with
		const uint8_t c,
		std::shared_ptr<FMIndex> pFM_index
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
	 * 
	 * @TODO: changed ik.start() < pFM_index->primary && ik.end() >= pFM_index->primary
	 * to current version... how to check if thats ok?
	 */
	if(
			ik.start() <= pFM_index->primary && 
			ik.end() > pFM_index->primary
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

    assert(cnts[c] >= 0);


	//BWAs SA intervals seem to be (a,b] while mine are [a,b)
	//pFM_index->L2[c] start of nuc c in BWT
	//cntk[c] offset of new interval
	//cntl[c] end of new interval
	return SAInterval(pFM_index->L2[c] + cntk[c] + 1, cntk_2[complement(c)], cnts[c]);
} // method


Interval<nucSeqIndex> BinarySeeding::lrExtension(
		nucSeqIndex center,
		std::shared_ptr<FMIndex> pFM_index,
		std::shared_ptr<NucSeq> pQuerySeq,
		std::shared_ptr<SegmentVector> pSegmentVector
	)
{
	// query sequence itself 
	const uint8_t *q = pQuerySeq->pGetSequenceRef(); 
	
	/* Initialize ik on the foundation of the single base q[x].
	 * In order to understand this initialization you should have a look 
	 * to the corresponding PowerPoint slide.
	 */
	// start I(q[x]) in T (start in BWT used for backward search) + 1, 
	// because very first string in SA-array starts with $
	// size in T and T' is equal due to symmetry
	SAInterval ik(
						pFM_index->L2[complement(q[center])] + 1, 
						pFM_index->L2[(int)q[center]] + 1, 
						pFM_index->L2[(int)q[center] + 1] - pFM_index->L2[(int)q[center]]
					);

	/*
	 * extend ik right, until there are no more matches
	 */
	nucSeqIndex end = center;
	for(nucSeqIndex i = center+1; i < pQuerySeq->length(); i++)
	{
		DEBUG_2(
			std::cout << i-1 << " -> " << ik.start() << " " << ik.end() << std::endl;
			std::cout << i-1 << " ~> " << ik.revComp().start() << " " << ik.revComp().end() << std::endl;
		)
		assert(ik.size() > 0);
		SAInterval ok = extend_backward(ik, complement(q[i]), pFM_index);

		DEBUG_2(
			std::cout << i << " -> " << ok.start() << " " << ok.end() << std::endl;
			std::cout << i << " ~> " << ok.revComp().start() << " " << ok.revComp().end() << std::endl;
		)
		/*
		* In fact, if ok.getSize is zero, then there are no matches any more.
		*/
		if (ok.size() == 0)
			break; // the SA-index interval size is too small to be extended further
		end = i;
		ik = ok;
	}//for
	DEBUG_2(
		std::cout << "swap" << std::endl;
	)
	//this is required in order to extend the other way
	ik = ik.revComp();
	nucSeqIndex start = center;
	/*
	 * extend ik left, until there are no more matches
	 */
	if(center > 0)
	{
		for(nucSeqIndex i = center-1; i >= 0; i--)
		{
			DEBUG_2(
				std::cout << i+1 << " -> " << ik.start() << " " << ik.end() << std::endl;
				std::cout << i+1 << " ~> " << ik.revComp().start() << " " << ik.revComp().end() << std::endl;
			)
			assert(ik.size() > 0);
			SAInterval ok = extend_backward(ik, q[i], pFM_index);
			DEBUG_2(
				std::cout << i << " -> " << ok.start() << " " << ok.end() << std::endl;
				std::cout << i << " ~> " << ok.revComp().start() << " " << ok.revComp().end() << std::endl;
			)

			/*
			* In fact, if ok.getSize is zero, then there are no matches any more.
			*/
			if (ok.size() == 0)
				break; // the SA-index interval size is too small to be extended further
			start = i;
			ik = ok;
			//cause nuxSeqIndex is unsigned
			if(i == 0)
				break;
		}//for
	}//if
	std::shared_ptr<Segment> pRightLeft(new Segment(start,end-start,ik));
	assert(start >= 0);
	assert(end < pQuerySeq->length());
	assert(pRightLeft->end() < pQuerySeq->length());
	pSegmentVector->push_back(pRightLeft);
	DEBUG_2(
		std::cout << "--other way--" << std::endl;
	)
	/* Initialize ik on the foundation of the single base q[x].
	 * In order to understand this initialization you should have a look 
	 *to the corresponding PowerPoint slide.
	 */
	// start I(q[x]) in T (start in BWT used for backward search) + 1, 
	// because very first string in SA-array starts with $
	// size in T and T' is equal due to symmetry
	ik = SAInterval(
						pFM_index->L2[q[center]] + 1, 
						pFM_index->L2[(int)complement(q[center])] + 1, 
						pFM_index->L2[(int)q[center] + 1] - pFM_index->L2[(int)q[center]]
					);
	start = center;
	/*
	 * extend ik left, until there are no more matches
	 */
	if(center > 0)
	{
		for(nucSeqIndex i = center-1; i >= 0; i--)
		{
			DEBUG_2(
				std::cout << i+1 << " -> " << ik.start() << " " << ik.end() << std::endl;
				std::cout << i+1 << " ~> " << ik.revComp().start() << " " << ik.revComp().end() << std::endl;
			)
			assert(ik.size() > 0);
			SAInterval ok = extend_backward(ik, q[i], pFM_index);
			DEBUG_2(
				std::cout << i << " -> " << ok.start() << " " << ok.end() << std::endl;
				std::cout << i << " ~> " << ok.revComp().start() << " " << ok.revComp().end() << std::endl;
			)

			/*
			* In fact, if ok.getSize is zero, then there are no matches any more.
			*/
			if (ok.size() == 0)
				break; // the SA-index interval size is too small to be extended further
			start = i;
			ik = ok;
			//cause nuxSeqIndex is unsigned
			if(i == 0)
				break;
		}//for
	}//if
	DEBUG_2(
		std::cout << "swap" << std::endl;
	)
	//this is required in order to extend the other way
	ik = ik.revComp();
	end = center;
	/*
	 * extend ik right, until there are no more matches
	 */
	for(nucSeqIndex i = center+1; i < pQuerySeq->length(); i++)
	{
		DEBUG_2(
			std::cout << i-1 << " -> " << ik.start() << " " << ik.end() << std::endl;
			std::cout << i-1 << " ~> " << ik.revComp().start() << " " << ik.revComp().end() << std::endl;
		)
		assert(ik.size() > 0);
		SAInterval ok = extend_backward(ik, complement(q[i]), pFM_index);

		DEBUG_2(
			std::cout << i << " -> " << ok.start() << " " << ok.end() << std::endl;
			std::cout << i << " ~> " << ok.revComp().start() << " " << ok.revComp().end() << std::endl;
		)

		/*
		* In fact, if ok.getSize is zero, then there are no matches any more.
		*/
		if (ok.size() == 0)
			break; // the SA-index interval size is too small to be extended further
		end = i;
		ik = ok;
	}//for
	std::shared_ptr<Segment> pLeftRight(new Segment(start,end-start,ik.revComp()));
	assert(start >= 0);
	assert(end < pQuerySeq->length());
	assert(pLeftRight->end() < pQuerySeq->length());
	pSegmentVector->push_back(pLeftRight);

	//to return the covered area
	Interval<nucSeqIndex> ret(center,0);
	if(pLeftRight->start() < pRightLeft->start())
		ret.start(pLeftRight->start());
	else
		ret.start(pRightLeft->start());

	if(pLeftRight->end() > pRightLeft->end())
		ret.end(pLeftRight->end());
	else
		ret.end(pRightLeft->end());

	return ret;
}//function

// @TODO: i should be my own module
// @TODO: test me
void BinarySeeding::bowtieExtension(
        std::shared_ptr<FMIndex> pFM_index,
        std::shared_ptr<NucSeq> pQuerySeq,
        std::shared_ptr<SegmentVector> pSegmentVector
    )
{
	const uint8_t *q = pQuerySeq->pGetSequenceRef(); 
    for(nucSeqIndex i = 0; i < pQuerySeq->length()-16; i+= 10)
    {

        SAInterval ik(
                            pFM_index->L2[complement(q[i])] + 1, 
                            pFM_index->L2[(int)q[i]] + 1, 
                            pFM_index->L2[(int)q[i] + 1] - pFM_index->L2[(int)q[i]]
                        );
        for(nucSeqIndex i2 = 0; i2 < 16; i2++)
        {
            // this is the extension
            // actually forward since we extend backwards on the reverse complement...
            ik = extend_backward(ik, complement(q[i2 + i]), pFM_index);
            if(ik.size() == 0)
                break;
        }//for
        if(ik.size() == 0)
            continue;
        // found a seed
        pSegmentVector->push_back(std::shared_ptr<Segment>(new Segment(i, 16, ik.revComp())));
    }//for
}//function

// @TODO: i should be my own module
// @TODO: test me 
// @TODO: figure out what the exactly do
/**
 * This is what blasr does:
 * 
 * for each position in the query:
 *  extend maximally backwards 
 *  save a seed that is one shorter than the maximal extension
 *  if the seed is longer than K = 12 adn maxambiguity = ?
 */
void BinarySeeding::doBlasrExtension(
        std::shared_ptr<FMIndex> pFM_index,
        std::shared_ptr<NucSeq> pQuerySeq,
        std::shared_ptr<SegmentVector> pSegmentVector
    )
{
	const uint8_t *q = pQuerySeq->pGetSequenceRef(); 
    for(nucSeqIndex i = 0; i < pQuerySeq->length(); i+= 1)
    {

        SAInterval ik(
                            pFM_index->L2[(int)q[i]] + 1, 
                            pFM_index->L2[complement(q[i])] + 1, 
                            pFM_index->L2[(int)q[i] + 1] - pFM_index->L2[(int)q[i]]
                        );
        SAInterval lk;// will hold the maximal extension
        SAInterval llk;// will hold the interval one shorter than the maximal extension
        nucSeqIndex i2 = 0;
        while(i2 <= i)
        {
            llk = lk;
            lk = ik;
            ik = extend_backward(ik, q[i - i2], pFM_index);
            if(ik.size() == 0)
                break;
            i2++;
        }//for
        if(i2 <= 12)
            continue;
        // found a seed
        pSegmentVector->push_back(std::shared_ptr<Segment>(new Segment(i - i2 + 1, i2-1, llk)));
    }//for
}//function

Interval<nucSeqIndex> BinarySeeding::nonEnclosedExtension(
		nucSeqIndex center,
		std::shared_ptr<FMIndex> pFM_index,
		std::shared_ptr<NucSeq> pQuerySeq,
		std::shared_ptr<SegmentVector> pSegmentVector
	)
{
	//to remember the covered area
	Interval<nucSeqIndex> ret(center,0);

	// query sequence itself
	const uint8_t *q = pQuerySeq->pGetSequenceRef(); 
	
	/* Initialize ik on the foundation of the single base q[x].
	 * In order to understand this initialization you should have a look 
	 *to the corresponding PowerPoint slide.
	 */
	// start I(q[x]) in T (start in BWT used for backward search) + 1, 
	// because very first string in SA-array starts with $
	// size in T and T' is equal due to symmetry
	SAInterval ik(
						pFM_index->L2[complement(q[center])] + 1, 
						pFM_index->L2[(int)q[center]] + 1, 
						pFM_index->L2[(int)q[center] + 1] - pFM_index->L2[(int)q[center]]
					);

	/*
	 * forward extension first
	 * this way we need to swap only once (forward to backwards) instead of swapping
	 * (backwards to forwards to backwards)
	 * 
	 * curr is used to remember the Suffix array interval each time we loose some hits by extending
	 */
	std::vector<Segment> curr;
	// extend until the end of the query
	for(nucSeqIndex i = center+1; i < pQuerySeq->length(); i++)
	{
		DEBUG_2(
			std::cout << i-1 << " -> " << ik.start() << " " << ik.end() << std::endl;
			std::cout << i-1 << " ~> " << ik.revComp().start() << " " << ik.revComp().end() << std::endl;
		)
		assert(ik.size() > 0);
		//this is the extension
		SAInterval ok = extend_backward(ik, complement(q[i]), pFM_index);

		// checking weather we lost some intervals
		// if so -> remember the interval just before we lost the hits
		if(ok.size() != ik.size()) 
			// save the reverse complement cause when extending the saved interval we will extend
			// in the other direction
			curr.push_back(Segment(center, i-center-1, ik.revComp()));
		// if were at the end of the query and we still have hits we need to make sure to record them
		if(i == pQuerySeq->length()-1 && ok.size() != 0)
			// save the reverse complement cause when extending the saved interval we will extend
			// in the other direction
			curr.push_back(Segment(center, i-center, ok.revComp()));

		DEBUG_2(
			std::cout << i << " -> " << ok.start() << " " << ok.end() << std::endl;
			std::cout << i << " ~> " << ok.revComp().start() << " " << ok.revComp().end() << std::endl;
		)
		/*
		* In fact, if ok.getSize is zero, then there are no matches any more.
		* thus we can stop extending forwards
		*/
		if (ok.size() == 0)
			break; // the SA-index interval size is too small to be extended further
		// if we get here we can forget the old interval and save the current interval.
		ik = ok;
		// remember that we covered this area
		ret.end(i);
	}//for
	DEBUG_2(
		std::cout << "swap" << std::endl;
	)
	/*
	 * This is the backwards extension part
	 * Here we need to extend the intervals in reverse order with respect to how we discovered them.
	 * (reversing is done by push_front insted of push_back)
	 *
	 * we will use prev and curr in this way:
	 * 		each iteration we will extend all intervals in prev
	 * 		and save the intervals that need to be extended further in curr
	 * 		at the end of the iteration we will swap prev and curr
	 * 		then clear curr
	 */
	std::reverse(curr.begin(), curr.end());
	std::vector<Segment> prev;
	//pointers for easy swapping of the lists

	// FIXME: for some reason valgrind does NOT like this
	// maybe it cant deal with the pointers?
	std::vector<Segment> *pPrev, *pCurr, *pTemp;
	pPrev = &curr;
	pCurr = &prev;
	// quick check that we can extend backwards at all (center is unsigned thus this is necessary)
	if(center != 0)
	{
		// extend until we reach the start of the query
		for(nucSeqIndex i = center-1; i >= 0; i--)
		{
			assert(pCurr->empty());

			/*
			 * we need to remember weather finished extending some interval in this step.
			 * because:
			 * 		if we already have found one with this length
			 * 			then all following intervals that we find have to be enclosed
			 * 			(this is due to the fact that they we know they start further right but 
			 * 			 end at the same point)
			 */
			bool bHaveOne = false;

			/*
			 * for all remembered intervals 
			 * (ordered by the start on the query)
			 */
			for(Segment& ik : *pPrev)
			{
				DEBUG_2(
					std::cout << i+1 << " -> " << ik.saInterval().start() << " " << ik.saInterval().end() << std::endl;
					std::cout << i+1 << " ~> " << ik.saInterval().revComp().start() << " " << ik.saInterval().revComp().end() << std::endl;
				)
				// actually extend the current interval
				SAInterval ok = extend_backward(ik.saInterval(), q[i], pFM_index);
				DEBUG_2(
					std::cout << i << " -> " << ok.start() << " " << ok.end() << std::endl;
					std::cout << i << " ~> " << ok.revComp().start() << " " << ok.revComp().end() << std::endl;
				)
				DEBUG(
					std::cout << ik.start() << ", " << ik.end() << ": " << ik.saInterval().size() << " -> " << ok.size() << std::endl;
				)
				// check if the extension resulted in a non enclosed interval
				if(ok.size() == 0 && !bHaveOne)
				{
					// save the interval
					pSegmentVector->push_back(std::shared_ptr<Segment>(new Segment(ik)));
					assert(ik.end() <= pQuerySeq->length());
					// we need to remember that we already found a interval this iteration
					bHaveOne = true;
				}// if
				// check if we can extend this interval further
				else if(ok.size() > 0)
				{
					// if so add the intervals to the list
					Segment xSeg = Segment(i, ik.size()+1, ok);
					// FIXME: memory leak here according to valgrind ?!?
					pCurr->push_back(xSeg);
					assert(xSeg.end() <= pQuerySeq->length());
				}// if
			}// for


			// swap out the lists and clear the things we just worked on
			pTemp = pPrev;
			pPrev = pCurr;
			pCurr = pTemp;
			// FIXME: memory leak here according to valgrind ?!?
			pCurr->clear();
			pCurr->shrink_to_fit();

			// if there are no more intervals to extend
			if(pPrev->empty())
				break;

			// remember that we covered this area
			ret.start(i);

			// cause nuxSeqIndex is unsigned we have to avoid underflow
			if(i == 0)
				break;
		}// for
	}// if

	//if we reach the beginning of the query it is possible that there are still intervals that contain matches.
	//we need to save the longest of those, which is conveniently (due to our sorting) the first one in the list
	if(!pPrev->empty())
	{
		assert(pPrev->front().size() >= pPrev->back().size());

		pSegmentVector->push_back(std::shared_ptr<Segment>(new Segment(pPrev->front())));
		assert(pPrev->front().end() <= pQuerySeq->length());
		
		DEBUG_2(
			std::cout << pPrev->front().start() << ":" << pPrev->front().end() << std::endl;
		)
	}//if

	//return the area that we covered
	return ret;
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
void BinarySeeding::procesInterval(
			Interval<nucSeqIndex> xAreaToCover,
			std::shared_ptr<SegmentVector> pSegmentVector,
			std::shared_ptr<FMIndex> pFM_index,
			std::shared_ptr<NucSeq> pQuerySeq,
			ThreadPoolAllowingRecursiveEnqueues* pxPool
		)
{
	nucSeqIndex uiStart = xAreaToCover.start();
	nucSeqIndex uiEnd = xAreaToCover.end();
	DEBUG(
		std::cout << "interval (" << uiStart << "," << uiEnd << ")" << std::endl;
	)

	Interval<nucSeqIndex> xAreaCovered;
	// performs extension and records any found seeds
	// here we use bLrExtension to choose the extension scheme
	if(bLrExtension)
		xAreaCovered = lrExtension(xAreaToCover.center(), pFM_index, pQuerySeq, pSegmentVector);
	else
		xAreaCovered = nonEnclosedExtension(xAreaToCover.center(), pFM_index, pQuerySeq, pSegmentVector);

	// extract how far the extension got on the query.
	nucSeqIndex uiFrom = xAreaCovered.start();
	nucSeqIndex uiTo = xAreaCovered.end();
	DEBUG(
		std::cout << "splitting interval (" << uiStart << "," << uiEnd << ") at (" << uiFrom << "," << uiTo << ")" << std::endl;
	)

	// if the extension did not fully cover until uiStart:
	if (uiFrom != 0 && uiStart + 1 < uiFrom)
	{
		// enqueue procesInterval() for a new interval that spans from uiStart to 
		// where the extension stopped
		pxPool->enqueue( 
			BinarySeeding::procesIntervalStatic,
			this,
			Interval<nucSeqIndex>(uiStart, uiFrom - uiStart - 1),
			pSegmentVector,
			pFM_index,
			pQuerySeq,
			pxPool
		);//enqueue
	}//if
	// if the extension did not fully cover until uiEnd:
	if (uiEnd > uiTo + 1)
	{
		// enqueue procesInterval() for a new interval that spans from where the extension stopped
		// to uiEnd
		pxPool->enqueue( 
			BinarySeeding::procesIntervalStatic,
			this,
			Interval<nucSeqIndex>(uiTo + 1, uiEnd - uiTo - 1),
			pSegmentVector,
			pFM_index,
			pQuerySeq,
			pxPool
		);//enqueue
	}//if
}//function


ContainerVector BinarySeeding::getInputType() const
{
	return ContainerVector{
			//the forward fm_index
			std::shared_ptr<Container>(new FMIndex()),
			//the query sequence
			std::shared_ptr<Container>(new NucSeq()),
		};
}
std::shared_ptr<Container> BinarySeeding::getOutputType() const
{
	return std::shared_ptr<Container>(new SegmentVector());
}


std::shared_ptr<Container> BinarySeeding::execute(
		std::shared_ptr<ContainerVector> vpInput
	)
{
	std::shared_ptr<FMIndex> pFM_index = std::static_pointer_cast<FMIndex>((*vpInput)[0]);
	std::shared_ptr<NucSeq> pQuerySeq = 
		std::static_pointer_cast<NucSeq>((*vpInput)[1]);

	std::shared_ptr<SegmentVector> pSegmentVector(new SegmentVector());

    if(do16ntevery10ntExtension)
		bowtieExtension(pFM_index, pQuerySeq, pSegmentVector);
    if(blasrExtension)
        doBlasrExtension(pFM_index, pQuerySeq, pSegmentVector);
    else
	{//scope for xPool
		// setup a threadpool
		ThreadPoolAllowingRecursiveEnqueues xPool( NUM_THREADS_ALIGNER );

		//enqueue the root interval (spanning the entire query) for processing
		xPool.enqueue( 
			BinarySeeding::procesIntervalStatic,
			this,
			Interval<nucSeqIndex>(0, pQuerySeq->length()),
			pSegmentVector,
			pFM_index,
			pQuerySeq,
			&xPool
		);//enqueue

	}//else & end of scope xPool

	return pSegmentVector;
}//function

void exportBinarySeeding()
{
	//export the BinarySeeding class
	boost::python::class_<
			BinarySeeding, 
			boost::python::bases<Module>,
        	std::shared_ptr<BinarySeeding>
		>(
			"BinarySeeding",
			boost::python::init<bool>()
		)
        .def_readwrite("do16ntevery10ntExtension", &BinarySeeding::do16ntevery10ntExtension)
        .def_readwrite("blasrExtension", &BinarySeeding::blasrExtension)
		;
	boost::python::implicitly_convertible< 
		std::shared_ptr<BinarySeeding>,
		std::shared_ptr<Module> 
	>();

}//function