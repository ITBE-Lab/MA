#include "segmentation.h"


std::shared_ptr<Container> Segmentation::execute(std::shared_ptr<Container> pInput)
{
    assert(*xSegmentTree.begin() != nullptr);

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
			xSegmentTree.begin(), this, &xPool
		);//enqueue

	}//end of scope xPool

	xSegmentTree.getTheNLongestIntervals(uiNumSegmentsAsAnchors).forEach(
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


#if confGENEREATE_ALIGNMENT_QUALITY_OUTPUT
	pxQuality->uiAmountSegments = xSegmentTree.length();
	pxQuality->uiLongestSegment = 0;
	nucSeqIndex uiEndLast = 0;
	if(*xSegmentTree.end() != nullptr)
		uiEndLast = xSegmentTree.end()->getStartIndex();
	pxQuality->uiLongestGapBetweenSegments = 0;
	xSegmentTree.forEach(
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

    return pInput;
}//function