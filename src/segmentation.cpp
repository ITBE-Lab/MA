#include "segmentation.h"

#define INCLUDE_SELF_TESTS (  1 )

//#include "BWTmem.h"
#include <vector>
//#include "analysis.h"
#include <memory>
#include <atomic>
#include <chrono>
//#include "assembly.h"
 

/* Emergency fix for broken forward search in original BWA code.
* Delivers 4 intervals for a single input interval.
* Here we use only two fields in the BWT_Interval.
* ik[0] == start I(P) in T	(start in BWT with backward search)
* ik[1] is unused
* ik[2] == size of interval (equal for I(P) and I'(P'), due to symmetry)
*/
SA_IndexInterval bwt_extend_backward(
									// current interval
									const SA_IndexInterval &ik,
									// the character to extend with
									const uint8_t c,
									// the FM index by reference
									std::shared_ptr<FM_Index> pxFM_Index
									)
{
	bwt64bitCounter cntk[4]; // Number of A, C, G, T in BWT until start of interval ik
	bwt64bitCounter cntl[4]; // Number of A, C, G, T in BWT until end of interval ik

	pxFM_Index->bwt_2occ4(
		// until start of SA index interval 
		// (-1, because we count the characters in front of the index)
		ik.start() - 1,
		// until end of SA index interval (-1, because bwt_occ4 counts inclusive)
		ik.end() - 1,
		cntk,						// output: Number of A, C, G, T until start of interval
		cntl						// output: Number of A, C, G, T until end of interval
	);
	//pxFM_Index->L2[c] start of nuc c in BWT
	//cntk[c] + 1 offset of new interval
	//cntl[c] length of new interval
	return SA_IndexInterval(pxFM_Index->L2[c] + 1 + cntk[c], cntl[c] - cntk[c]);
} // method

bool Segmentation::canExtendFurther(std::shared_ptr<SegmentTreeInterval> pxNode, nucSeqIndex uiCurrIndex, bool bBackwards, nucSeqIndex uiQueryLength)
{
	if (!bBackwards)
	{
		return uiCurrIndex < uiQueryLength;
	}//if
	else
	{
		return uiCurrIndex > 0;
	}//else

	return true;
	//we want to allow extension past the interval borders
	if (!bBackwards)
	{
		return uiCurrIndex <= pxNode->end();
	}//if
	else
	{
		return uiCurrIndex >= pxNode->start();
	}//else
}//function

bool isMinIntervalSizeReached(nucSeqIndex uiStartIndex, nucSeqIndex uiCurrIndex, 
							  nucSeqIndex uiMinIntervalSize)
{
	if (uiStartIndex > uiCurrIndex)
		return uiStartIndex - uiCurrIndex >= uiMinIntervalSize;
	else if (uiStartIndex < uiCurrIndex)
		return uiCurrIndex - uiStartIndex >= uiMinIntervalSize;
	return false;
}//function

bool isFurtherThan(
					bool bBackwards, nucSeqIndex uiCurrIndex, 
					nucSeqIndex uiOnlyRecordHitsFurtherThan
				  )
{
	if (bBackwards)
		return uiCurrIndex < uiOnlyRecordHitsFurtherThan;
	else
		return uiCurrIndex > uiOnlyRecordHitsFurtherThan;
}//function

nucSeqIndex Segmentation::extend(
								 std::shared_ptr<SegmentTreeInterval> pxNode, 
								 nucSeqIndex uiStartIndex, bool bBackwards, 
								 nucSeqIndex uiOnlyRecordHitsFurtherThan
								)
{
	nucSeqIndex i = uiStartIndex;
	assert(i >= 0 && i < pxQuerySeq->length());

	const uint8_t *q = pxQuerySeq->pGetSequenceRef(); // query sequence itself

	//select the correct fm_index based on the direction of the extension
	std::shared_ptr<FM_Index> pxUsedFmIndex;
	if (!bBackwards)
		pxUsedFmIndex = pxRev_FM_Index;
	else
		pxUsedFmIndex = pxFM_index;
	
	/* Initialize ik on the foundation of the single base q[x].
	 * In order to understand this initialization you should have a look 
	 *to the corresponding PowerPoint slide.
	 */
	// start I(q[x]) in T (start in BWT used for backward search) + 1, 
	// because very first string in SA-array starts with $
	// size in T and T' is equal due to symmetry
	SA_IndexInterval ik(
						pxUsedFmIndex->L2[(int)q[i]] + 1, 
						pxUsedFmIndex->L2[(int)q[i] + 1] - pxUsedFmIndex->L2[(int)q[i]]
					);

	if (!bBackwards)
		i++;
	else if (i == 0)// unsigned int
		return 0;
	else
		i--;

	// while we can extend further
	while ( canExtendFurther(pxNode ,i ,bBackwards ,pxQuerySeq->length()) )
	{
		assert(i >= 0 && i < pxQuerySeq->length());
		if (q[i] >= 4 && bBreakOnAmbiguousBase) // An ambiguous base
		{
			break; // break if the parameter is set.
		}//if

		const uint8_t c = q[i]; // character at position i in the query

		// perform one step of the extension
		SA_IndexInterval ok = bwt_extend_backward(ik, c, pxUsedFmIndex);

		/* 
		 * Pick the extension interval for character c 
		 * and check with respect to changing interval size.
		 */
		if (ok.size() != ik.size())
		{
			/*
			* In fact, if ok.getSize is zero, then there are no matches any more.
			*/
			if (ok.size() == 0)
			{
				break; // the SA-index interval size is too small to be extended further
			} // if

			// once the min interval size is reached, 
			// record the matches every time we lose some by extending further.
			if (
					isMinIntervalSizeReached(uiStartIndex, i, uiMinIntervalSize) && 
					isFurtherThan(bBackwards, i, uiOnlyRecordHitsFurtherThan)
				)
			{
				// record the current match
				if (!bBackwards)
					pxNode->push_back(SaSegment(uiStartIndex, i - uiStartIndex - 1, ik, true), false);
				else
					pxNode->push_back(SaSegment(i + 1, uiStartIndex - i - 1, ik, false), false);
			}//if
		}//if

		/* Set input interval for the next iteration.
		*/
		ik = ok;

		if (!bBackwards)
			i++;
		else if (i == 0)// unsigned int
			break;
		else
			i--;
	}//while


	//once the min interval size is reached, record the matches every time we lose some by extending further.
	if (
			isMinIntervalSizeReached(uiStartIndex, i, uiMinIntervalSize) && 
			isFurtherThan(bBackwards, i, uiOnlyRecordHitsFurtherThan)
		)
	{
		// record the current match
		if (!bBackwards)
			pxNode->push_back(SaSegment(uiStartIndex, i - uiStartIndex - 1, ik, true), true);
		else
			pxNode->push_back(SaSegment(i + 1, uiStartIndex - i - 1, ik, false), true);
	}//if

	if (!bBackwards)
		return i-1;
	else
		return i+1;
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
void Segmentation::procesInterval(size_t uiThreadId, SegTreeItt pxNode, ThreadPoolAllowingRecursiveEnqueues *pxPool)
{
	//performs backwards extension and records any perfect matches
	nucSeqIndex uiBackwardsExtensionReachedUntil = extend(*pxNode, pxNode->getCenter(), true, pxNode->getCenter());
	/* we could already extend our matches from the center of the node to  uiBackwardsExtensionReachedUntil
	*  therefore we know that the next extension will reach at least as far as the center of the node.
	*  we might be able to extend our matches further though.
	*/
	nucSeqIndex uiBackwardsForwardsExtensionReachedUntil = extend(*pxNode, uiBackwardsExtensionReachedUntil, false, uiBackwardsExtensionReachedUntil);
	assert(uiBackwardsForwardsExtensionReachedUntil >= pxNode->getCenter());

	//performs forward extension and records any perfect matches
	nucSeqIndex uiForwardsExtensionReachedUntil = extend(*pxNode, pxNode->getCenter(), false, pxNode->getCenter());
	/* we could already extend our matches from the center of the node to  uiForwardExtensionReachedUntil
	*  therefore we know that the next extension will reach at least as far as the center of the node.
	*  we might be able to extend our matches further though.
	*/
	nucSeqIndex uiForwardsBackwardsExtensionReachedUntil = extend(*pxNode, uiForwardsExtensionReachedUntil, true, uiForwardsExtensionReachedUntil);
	assert(uiForwardsBackwardsExtensionReachedUntil <= pxNode->getCenter());

	nucSeqIndex uiFrom, uiTo;
	if (uiBackwardsForwardsExtensionReachedUntil - uiBackwardsExtensionReachedUntil > uiForwardsExtensionReachedUntil - uiForwardsBackwardsExtensionReachedUntil)
	{
		uiFrom = uiBackwardsExtensionReachedUntil;
		uiTo = uiBackwardsForwardsExtensionReachedUntil;
	}//if
	else
	{
		uiFrom = uiForwardsBackwardsExtensionReachedUntil;
		uiTo = uiForwardsExtensionReachedUntil;
	}//else

	//if the prev interval is longer than uiMinIntervalSize
	if (pxNode->start() + uiMinIntervalSize <= uiFrom)
	{
		//create a new list element and insert it before the current node
		auto pxPrevNode = pSegmentTree->insertBefore(std::shared_ptr<SegmentTreeInterval>(
			new SegmentTreeInterval(pxNode->start(), uiFrom - pxNode->start())), pxNode);
		//enqueue procesInterval() for the new interval
		pxPool->enqueue(
			[/*WARNING: do not catch anything here: the lambda function is enqueued to into a thread pool, 
			 local variables might not exist anymore during it's execution*/]
			(size_t uiThreadId, SegTreeItt pxPrevNode, Segmentation *pxAligner, ThreadPoolAllowingRecursiveEnqueues *pxPool)
			{
				pxAligner->procesInterval(uiThreadId, pxPrevNode, pxPool);
			},//lambda
			pxPrevNode, this, pxPool
		);//enqueue
	}//if
	//if the post interval is longer than uiMinIntervalSize
	if (pxNode->end() >= uiTo + uiMinIntervalSize)
	{
		//create a new list element and insert it after the current node
		auto pxNextNode = pSegmentTree->insertAfter(std::shared_ptr<SegmentTreeInterval>(
			new SegmentTreeInterval(uiTo, pxNode->end() - uiTo)), pxNode);
		//enqueue procesInterval() for the new interval
		pxPool->enqueue(
			[/*WARNING: do not catch anything here: the lambda function is enqueued to into a thread pool,
			 local variables might not exist anymore during it's execution*/]
			(size_t uiThreadId, SegTreeItt pxNextNode, Segmentation *pxAligner, ThreadPoolAllowingRecursiveEnqueues *pxPool)
			{
				pxAligner->procesInterval(uiThreadId, pxNextNode, pxPool);
			},//lambda
			pxNextNode, this, pxPool
		);//enqueue
	}//if


	DEBUG(
		std::cout << "splitting interval " << *pxNode << " at (" << uiFrom << "," << uiTo << ")" << std::endl;
	)
	pxNode->start(uiFrom);
	pxNode->end(uiTo);

	if (pxNode->size() < uiMinIntervalSize)
		pSegmentTree->removeNode(pxNode);
}//function

void Segmentation::segment()
{
	assert(*pSegmentTree->begin() != nullptr);

	{//scope for xPool
		ThreadPoolAllowingRecursiveEnqueues xPool(NUM_THREADS_ALIGNER);

		//enqueue the root interval for processing
		xPool.enqueue(
			[/*WARNING: do not catch anything here: the lambda function is enqueued to into a thread pool,
			 local variables might not exist anymore during it's execution*/]
			(size_t uiThreadId, SegTreeItt pxRoot, Segmentation *pxAligner, ThreadPoolAllowingRecursiveEnqueues *pxPool)
			{
				pxAligner->procesInterval(uiThreadId, pxRoot, pxPool);
			},//lambda
			pSegmentTree->begin(), this, &xPool
		);//enqueue

	}//end of scope xPool
}//function


std::vector<ContainerType> SegmentationContainer::getInputType()
{
	return std::vector<ContainerType>{
			//the forward fm_index
			ContainerType::fM_index,
			//the reversed fm_index
			ContainerType::fM_index,
			//the querry sequence
			ContainerType::nucSeq,
			//the reference sequence (packed since it could be really long)
			ContainerType::packedNucSeq,
		};
}
std::vector<ContainerType> SegmentationContainer::getOutputType()
{
	return std::vector<ContainerType>{ContainerType::segmentList};
}


std::shared_ptr<Container> SegmentationContainer::execute(
		std::vector<std::shared_ptr<Container>> vpInput
	)
{
	std::shared_ptr<FM_Index> pFM_index = std::static_pointer_cast<FM_Index>(vpInput[0]);
	std::shared_ptr<FM_Index> pFM_indexReversed = std::static_pointer_cast<FM_Index>(vpInput[1]);
	std::shared_ptr<NucleotideSequence> pQuerrySeq = 
		std::static_pointer_cast<NucleotideSequence>(vpInput[2]);
	std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pRefSeq = 
		std::static_pointer_cast<BWACompatiblePackedNucleotideSequencesCollection>(vpInput[3]);


	Segmentation xS(
			pFM_index, 
			pFM_indexReversed, 
			pQuerrySeq, 
			bBreakOnAmbiguousBase,
			bSkipLongBWTIntervals,
			uiMinIntervalSize,
			uiMaxHitsPerInterval,
			pRefSeq
		);
	xS.segment();

	return xS.pSegmentTree;
}//function


void exportSegmentation()
{
	//export the segmentation class
	boost::python::class_<SegmentationContainer, boost::python::bases<Module>>(
			"Segmentation",
			"bBreakOnAmbiguousBase: weather the extension of "
			"intervalls shall be stopped at N's\n"
			"bSkipLongBWTIntervals: skip seeds that have more than "
			"uiMaxHitsPerInterval matches on the reference\n"
			"uiMinIntervalSize: only record intervals greater than uiMinIntervalSize\n"
			"uiMaxHitsPerInterval: skip seeds that have more than "
			"uiMaxHitsPerInterval matches on the reference if bSkipLongBWTIntervals is set\n",
			boost::python::init<boost::python::optional<bool, bool, nucSeqIndex, unsigned int>>(
				"arg1: self\n"
				"arg2: weather the extension of "
				"intervalls shall be stopped at N's\n"
				"arg3: skip seeds that have more than "
				"arg5 matches on the reference\n"
				"arg4: only record intervals greater than arg4\n"
				"arg5: skip seeds that have more than "
				"arg5 matches on the reference if arg3 is set\n"
			)
		)
		.def_readwrite("bBreakOnAmbiguousBase", &SegmentationContainer::bBreakOnAmbiguousBase)
		.def_readwrite("bSkipLongBWTIntervals", &SegmentationContainer::bSkipLongBWTIntervals)
		.def_readwrite("uiMinIntervalSize", &SegmentationContainer::uiMinIntervalSize)
		.def_readwrite("uiMaxHitsPerInterval", &SegmentationContainer::uiMaxHitsPerInterval)
		;
	
}//function