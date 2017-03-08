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
void bwt_extend_backward(const SA_IndexInterval &ik,	// (input) single interval
						SA_IndexInterval ok[4],		// (output) 4 intervals, for each symbol 1 (ok[0] <-> A, ok[1] <-> C, ok[2] <-> G, ok[3] <-> T
						std::shared_ptr<FM_Index> pxFM_Index // the FM index by reference
					)
{
	//// xSavePrint << "bwt_extend with ik: " << ik.endIndexOfMatchInQuery() << " [" << ik.x[0] << " ," << ik.x[1] << " ," << ik.x[2] << "]\n";
	bwt64bitCounter cntk[4]; // Number of A, C, G, T in BWT until start of interval ik
	bwt64bitCounter cntl[4]; // Number of A, C, G, T in BWT until end of interval ik

	pxFM_Index->bwt_2occ4(ik.x[0] - 1,				// until start of SA index interval (-1, because we count the characters in front of the index)
		ik.x[0] + ik.x[2] - 1,	// until end of SA index interval (-1, because bwt_occ4 counts inclusive)
		cntk,						// output: Number of A, C, G, T until start of interval
		cntl						// output: Number of A, C, G, T until end of interval
		);

	/* Standard backward extension for the 4 symbols A, C, G, T.
	* So, we get I(cP) in ok[c] for all c in {A, C, G, T}.
	*/
	for (int c = 0; c < 4; ++c)
	{
		//// xSavePrint << "xFM_Index.L2[c] is " << c << " " << xFM_Index.L2[c] << " " << cntk[c] << "" << cntl[c] << "\n";
		ok[c].x[0] = pxFM_Index->L2[c] + 1 + cntk[c]; // start of Interval I(cP)  

		ok[c].x[2] = cntl[c] - cntk[c]; // size of Interval I(cP) == size of Interval I'(P'c')
		//// xSavePrint << "ok[c].x[2] is " << ok[c].x[2] << " " << cntl[c] << " " << cntk[c] << "\n";
	} // for
} // method

bool Segmentation::canExtendFurther(std::shared_ptr<SegmentTreeInterval> pxNode, nucSeqIndex uiCurrIndex, bool bBackwards, nucSeqIndex uiQueryLength)
{
	if (!bBackwards)
	{
		return uiCurrIndex < uiQueryLength;
	}//if
	else
	{
		return uiCurrIndex >= 0;
	}//else

	return true;
	//we want to allow extension past the interval borders
	if (!bBackwards)
	{
		return uiCurrIndex <= pxNode->getEndIndex();
	}//if
	else
	{
		return uiCurrIndex >= pxNode->getStartIndex();
	}//else
}//function

bool isMinIntervalSizeReached(nucSeqIndex uiStartIndex, nucSeqIndex uiCurrIndex, nucSeqIndex uiMinIntervalSize)
{
	if (uiStartIndex > uiCurrIndex)
		return uiStartIndex - uiCurrIndex >= uiMinIntervalSize;
	else if (uiStartIndex < uiCurrIndex)
		return uiCurrIndex - uiStartIndex >= uiMinIntervalSize;
	return false;
}//function

bool isFurtherThan(bool bBackwards, nucSeqIndex uiCurrIndex, nucSeqIndex uiOnlyRecordHitsFurtherThan)
{
	if (bBackwards)
		return uiCurrIndex < uiOnlyRecordHitsFurtherThan;
	else
		return uiCurrIndex > uiOnlyRecordHitsFurtherThan;
}//function

nucSeqIndex Segmentation::extend(std::shared_ptr<SegmentTreeInterval> pxNode, nucSeqIndex uiStartIndex, bool bBackwards, nucSeqIndex uiOnlyRecordHitsFurtherThan)
{
	nucSeqIndex i = uiStartIndex;
	assert(i >= 0 && i < pxQuerySeq->length());

	const uint8_t *q = pxQuerySeq->pGetSequenceRef(); // query sequence itself

	/* BIDI BWT intervals. Input interval and 4 output intervals for A/C/G/T.
	*/
	SA_IndexInterval ik, ok[4];

	//select the correct fm_index based on the direction of the extension
	std::shared_ptr<FM_Index> pxUsedFmIndex;
	if (!bBackwards)
		pxUsedFmIndex = pxRev_FM_Index;
	else
		pxUsedFmIndex = pxFM_index;
	
	/* Initialize ik on the foundation of the single base q[x].
	* Original code: bwt_set_intv(bwt, q[x], ik);
	* In order to understand this initialization you should have a look to the corresponding PowerPoint slide.
	*/
	ik.x[0] = pxUsedFmIndex->L2[(int)q[i]] + 1;								// start I(q[x]) in T (start in BWT used for backward search) + 1, because very first string in SA-array starts with $
	ik.x[2] = pxUsedFmIndex->L2[(int)q[i] + 1] - pxUsedFmIndex->L2[(int)q[i]];	// size in T and T' is equal due to symmetry

	if (!bBackwards)
		i++;
	else if (i == 0)// unsigned int
		return 0;
	else
		i--;

	//std::cout << "extending: (" << i << ") " << std::endl;

	//while we can extend further
	while (canExtendFurther(pxNode,i,bBackwards,pxQuerySeq->length()))
	{
		assert(i >= 0 && i < pxQuerySeq->length());
		if (q[i] >= 4 && bBreakOnAmbiguousBase) // An ambiguous base
		{
#ifdef DEBUG_ALIGNER
			BOOST_LOG_TRIVIAL( info ) << "break on ambiguous base";
#endif
			break; // break if the parameter is set.
		}//if

		//perform one step of the extension
		bwt_extend_backward(ik, ok, pxUsedFmIndex);

		const uint8_t c = q[i]; // character at position i in the query

#if 0
		NucleotideSequence xSec;
		auto iBegin = pxUsedFmIndex->bwt_sa(ik.x[0]);
		auto iEnd = iBegin + (bBackwards ? uiStartIndex - i : i - uiStartIndex);
		std::cout << "query: " << *pxQuerySeq->fromTo(bBackwards ? i + 1 : uiStartIndex, !bBackwards ? i : uiStartIndex + 1).asSequenceOfACGT() << " " << ik.x[2] << std::endl;
		if (bBackwards)
			pxRefSequence->vExtractSubsection(iBegin, iEnd, xSec);
		else
			pxRefSequence->vExtractSubsection(
					pxRefSequence->uiUnpackedSizeForwardPlusReverse() - iEnd, 
					pxRefSequence->uiUnpackedSizeForwardPlusReverse() - iBegin,
					xSec
				);
		
		std::cout << "reference:" << *xSec.fullSequenceAsSlice().asSequenceOfACGT() << std::endl;
		std::cout << iBegin << " " << iEnd << " \n"
			<< pxRefSequence->uiUnpackedSizeForwardPlusReverse() - iEnd << " " << pxRefSequence->uiUnpackedSizeForwardPlusReverse() - iBegin << std::endl;
#endif
		

		/* Pick the extension interval for character c and check with respect to changing interval size.
		*/
		if (ok[c].x[2] != ik.x[2])
		{
			/*
			* In fact, if ok[c].x[2] is zero, then there are no matches any more.
			*/
			if (ok[c].x[2] == 0)
			{
				break; // the SA-index interval size is too small to be extended further
			} // if

			//once the min interval size is reached, record the matches every time we lose some by extending further.
			if (isMinIntervalSizeReached(uiStartIndex, i, uiMinIntervalSize) && isFurtherThan(bBackwards, i, uiOnlyRecordHitsFurtherThan))
			{
				//record the current match
				if (!bBackwards)
					pxNode->pushBackBwtInterval(ik.x[0], ik.x[2], uiStartIndex, i, true, false);
				else
					pxNode->pushBackBwtInterval(ik.x[0], ik.x[2], i + 1, uiStartIndex + 1, false, false);
			}//if
		}//if

		/* Set input interval for the next iteration.
		*/
		ik = ok[c];

		if (!bBackwards)
			i++;
		else if (i == 0)// unsigned int
			break;
		else
			i--;
	}//while


	//once the min interval size is reached, record the matches every time we lose some by extending further.
	if (isMinIntervalSizeReached(uiStartIndex, i, uiMinIntervalSize) && isFurtherThan(bBackwards, i, uiOnlyRecordHitsFurtherThan))
	{
		//record the current match
		if (!bBackwards)
			pxNode->pushBackBwtInterval(ik.x[0], ik.x[2], uiStartIndex, i, true, true);
		else
			pxNode->pushBackBwtInterval(ik.x[0], ik.x[2], i + 1, uiStartIndex + 1, false, true);
	}//if

	//std::cout << " (" << i + 1 << ")" << std::endl;
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
	if (pxNode->getStartIndex() + uiMinIntervalSize <= uiFrom)
	{
		//create a new list element and insert it before the current node
		auto pxPrevNode = pSegmentTree->insertBefore(std::shared_ptr<SegmentTreeInterval>(
			new SegmentTreeInterval(pxNode->getStartIndex(), uiFrom)), pxNode);
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
	if (pxNode->getEndIndex() >= uiTo + uiMinIntervalSize)
	{
		//create a new list element and insert it after the current node
		auto pxNextNode = pSegmentTree->insertAfter(std::shared_ptr<SegmentTreeInterval>(
			new SegmentTreeInterval(uiTo, pxNode->getEndIndex())), pxNode);
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
#ifdef DEBUG_ALIGNER
	BOOST_LOG_TRIVIAL(info) << "splitting interval " << *pxNode << " at (" << uiFrom << "," << uiTo << ")";
#endif //DEBUG_ALIGNER
	pxNode->setInterval(uiFrom, uiTo);

	/*
	 *	try to perform backwards extension from 3/4 of the interval
	 *	as well as forwards extension from 1/4 of the interval
	 */
	nucSeqIndex uiThreeQuaterOfTheInterval = pxNode->getStartIndex() + pxNode->length() * 3 / 4;
	nucSeqIndex uiOneQuaterOfTheInterval = pxNode->getStartIndex() + pxNode->length() / 4;

	nucSeqIndex uiReachedUntil = extend(*pxNode, uiOneQuaterOfTheInterval, false, uiOneQuaterOfTheInterval);
	assert(uiReachedUntil >= pxNode->getCenter());

	uiReachedUntil = extend(*pxNode, uiThreeQuaterOfTheInterval, true, uiThreeQuaterOfTheInterval);
	assert(uiReachedUntil <= pxNode->getCenter());


	//if (pxNode->length() > uiMinIntervalSize)
	//	saveHits(*pxNode, uiThreadId); deprecated see @saveHits
	//else
	if (pxNode->length() < uiMinIntervalSize)
		pSegmentTree->removeNode(pxNode);

	//TODO: if we need to generate more hits we could split intervals longer than x even if we find matches in them, just to search from different starting positions.
	//in this case we should limit the splitting to n numbers of iterations though, otherwise the procedure will degenerate to a bwa-like search
}//function


void Segmentation::forEachNonBridgingHitOnTheRefSeq(std::shared_ptr<SegmentTreeInterval> pxNode, bool bAnchorOnly,
	std::function<void(nucSeqIndex ulIndexOnRefSeq, nucSeqIndex uiQueryBegin, nucSeqIndex uiQueryEnd)> fDo)
{
	pxNode->forEachHitOnTheRefSeq(
		pxFM_index, pxRev_FM_Index, uiMaxHitsPerInterval, bSkipLongBWTIntervals, bAnchorOnly,
#if confGENEREATE_ALIGNMENT_QUALITY_OUTPUT
	    pxQuality,
#endif
		[&](nucSeqIndex ulIndexOnRefSeq, nucSeqIndex uiQuerryBegin, nucSeqIndex uiQuerryEnd)
		{
			int64_t iSequenceId;
			//check if the match is bridging the forward/reverse strand or bridging between two chromosomes
			/* we have to make sure that the match does not start before or end after the reference sequence
			* this can happen since we can find parts on the end of the query at the very beginning of the reference or vis versa.
			* in this case we will replace the out of bounds index with 0 or the length of the reference sequence respectively.
			*/
			if (pxRefSequence->bridingSubsection(
				ulIndexOnRefSeq > uiQuerryBegin ? (uint64_t)ulIndexOnRefSeq - (uint64_t)uiQuerryBegin : 0,
				ulIndexOnRefSeq + pxQuerySeq->length() >= pxFM_index->getRefSeqLength() + uiQuerryBegin ? pxFM_index->getRefSeqLength() - ulIndexOnRefSeq : pxQuerySeq->length(),
				iSequenceId)
				)
			{
#ifdef DEBUG_CHECK_INTERVALS
			BOOST_LOG_TRIVIAL(info) << "skipping hit on bridging section (" << ulIndexOnRefSeq - uiQuerryBegin << ") for the interval " << *pxNode;
#endif
				//if so ignore this hit
				return;
			}//if
			fDo(ulIndexOnRefSeq, uiQuerryBegin, uiQuerryEnd);
		}//lambda
	);
}//function

void Segmentation::forEachNonBridgingPerfectMatch(std::shared_ptr<SegmentTreeInterval> pxNode, bool bAnchorOnly,
	std::function<void(std::shared_ptr<PerfectMatch>)> fDo)
{
	forEachNonBridgingHitOnTheRefSeq(
		pxNode, bAnchorOnly,
		[&](nucSeqIndex ulIndexOnRefSeq, nucSeqIndex uiQuerryBegin, nucSeqIndex uiQuerryEnd)
		{
			fDo(std::shared_ptr<PerfectMatch>(new PerfectMatch(uiQuerryEnd - uiQuerryBegin, ulIndexOnRefSeq, uiQuerryBegin)));
		}//lambda
	);//for each
}//function

#if 0
/*
*	deprecated since the anchor matches will be extracted after the segmentation process.
* 	the finding achor matches process will have it's own module
*/

/* transfer the saved hits into the clustering
 * if DEBUG_CHECK_INTERVALS is activated the hits are verified before storing
*/
void Segmentation::saveHits(std::shared_ptr<SegmentTreeInterval> pxNode, size_t uiThreadId)
{
	forEachNonBridgingPerfectMatch(
		pxNode, false,
		[&](std::shared_ptr<PerfectMatch> pxMatch)
		{
			pxAnchorMatchList->addMatch(std::shared_ptr<PerfectMatch>(pxMatch), uiThreadId);
		}//lambda
	);//for each
}///function

#endif

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

	/* TODO: move me into my own module
	pSegmentTree->getTheNLongestIntervals(uiNumSegmentsAsAnchors).forEach(
		[&](std::shared_ptr<SegmentTreeInterval> pxInterval)
		{
			forEachNonBridgingPerfectMatch(pxInterval, true,
				[&](std::shared_ptr<PerfectMatch> pxMatch)
				{
					pxAnchorMatchList->addAnchorSegment(pxMatch);
				}//lambda
			);
		}//lambda
	);
	*/

#if confGENEREATE_ALIGNMENT_QUALITY_OUTPUT
	pxQuality->uiAmountSegments = pSegmentTree->length();
	pxQuality->uiLongestSegment = 0;
	nucSeqIndex uiEndLast = 0;
	if(*pSegmentTree->end() != nullptr)
		uiEndLast = pSegmentTree->end()->getStartIndex();
	pxQuality->uiLongestGapBetweenSegments = 0;
	pSegmentTree->forEach(
		[&](std::shared_ptr<SegmentTreeInterval> pxInterval)
		{
			if (pxQuality->uiLongestGapBetweenSegments < pxInterval->getStartIndex() - uiEndLast)
				pxQuality->uiLongestGapBetweenSegments = pxInterval->getStartIndex() - uiEndLast;
			if (pxInterval->length() > pxQuality->uiLongestSegment)
				pxQuality->uiLongestSegment = pxInterval->length();
			uiEndLast = pxInterval->getEndIndex();
		}//lambda
	);
	pxQuality->uiAmountAnchorSegments = uiNumSegmentsAsAnchors;
#endif
}//function


std::shared_ptr<Container> SegmentationContainer::getInputType()
{
	std::shared_ptr<ContainerVector> pRet(new ContainerVector());
	//the forward fm_index
	pRet->vElements.push_back(std::shared_ptr<Container>(new DummyContainer(ContainerType::fM_index)));
	//the reversed fm_index
	pRet->vElements.push_back(std::shared_ptr<Container>(new DummyContainer(ContainerType::fM_index)));
	//the querry sequence
	pRet->vElements.push_back(std::shared_ptr<Container>(new DummyContainer(ContainerType::nucSeq)));
	//the reference sequence (packed since it could be really long)
	pRet->vElements.push_back(std::shared_ptr<Container>(new DummyContainer(ContainerType::packedNucSeq)));
	return pRet;
}
std::shared_ptr<Container> SegmentationContainer::getOutputType()
{
	return std::shared_ptr<Container>(new DummyContainer(ContainerType::segmentList));
}


std::shared_ptr<Container> SegmentationContainer::execute(std::shared_ptr<Container> pInput)
{

	std::shared_ptr<ContainerVector> pCastedInput = std::static_pointer_cast<ContainerVector>(pInput);
	std::shared_ptr<FM_IndexContainer> pFM_index = std::static_pointer_cast<FM_IndexContainer>(pCastedInput->vElements.at(0));
	std::shared_ptr<FM_IndexContainer> pFM_indexReversed = std::static_pointer_cast<FM_IndexContainer>(pCastedInput->vElements.at(1));
	std::shared_ptr<NucSeqContainer> pQuerrySeq = std::static_pointer_cast<NucSeqContainer>(pCastedInput->vElements.at(2));
	std::shared_ptr<PackContainer> pRefSeq = std::static_pointer_cast<PackContainer>(pCastedInput->vElements.at(3));


	Segmentation xS(pFM_index->pIndex, pFM_indexReversed->pIndex, pQuerrySeq->pSeq, bBreakOnAmbiguousBase, bSkipLongBWTIntervals, uiMinIntervalSize, uiMaxHitsPerInterval, pRefSeq->pPack);
	xS.segment();

	return std::shared_ptr<SegmentTreeContainer>(new SegmentTreeContainer(xS.pSegmentTree));
}//function


void exportSegmentation()
{

    //export the segmentation class
	boost::python::class_<SegmentationContainer, boost::python::bases<Module>>("Segmentation",boost::python::init<boost::python::optional<bool, bool, nucSeqIndex, unsigned int>>())
		.add_property("bBreakOnAmbiguousBase", &SegmentationContainer::bBreakOnAmbiguousBase, &SegmentationContainer::bBreakOnAmbiguousBase)
		.add_property("bSkipLongBWTIntervals", &SegmentationContainer::bSkipLongBWTIntervals, &SegmentationContainer::bSkipLongBWTIntervals)
		.add_property("uiMinIntervalSize", &SegmentationContainer::uiMinIntervalSize, &SegmentationContainer::uiMinIntervalSize)
		.add_property("uiMaxHitsPerInterval", &SegmentationContainer::uiMaxHitsPerInterval, &SegmentationContainer::uiMaxHitsPerInterval)
		;

}//function