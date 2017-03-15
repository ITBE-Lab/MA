#include "graphicalMethod.h"


std::shared_ptr<Container> Bucketing::getInputType()
{
	std::shared_ptr<ContainerVector> pRet(new ContainerVector());
	//all segments
	pRet->vElements.push_back(std::shared_ptr<Container>(new DummyContainer(ContainerType::segmentList)));
	//the anchors
	pRet->vElements.push_back(std::shared_ptr<Container>(new DummyContainer(ContainerType::segmentList)));
	//the querry
	pRet->vElements.push_back(std::shared_ptr<Container>(new DummyContainer(ContainerType::nucSeq)));
	//the reference
	pRet->vElements.push_back(std::shared_ptr<Container>(new DummyContainer(ContainerType::packedNucSeq)));
	//the forward fm_index
	pRet->vElements.push_back(std::shared_ptr<Container>(new DummyContainer(ContainerType::fM_index)));
	//the reversed fm_index
	pRet->vElements.push_back(std::shared_ptr<Container>(new DummyContainer(ContainerType::fM_index)));
	return pRet;
}//function

std::shared_ptr<Container> Bucketing::getOutputType()
{
	return std::shared_ptr<Container>(new DummyContainer(ContainerType::stripOfConsiderationList));
}//function


void Bucketing::forEachNonBridgingHitOnTheRefSeq(std::shared_ptr<SegmentTreeInterval> pxNode, bool bAnchorOnly, std::shared_ptr<FM_Index> pxFM_index, std::shared_ptr<FM_Index> pxRev_FM_Index, std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence, std::shared_ptr<NucleotideSequence> pxQuerySeq,
	std::function<void(nucSeqIndex ulIndexOnRefSeq, nucSeqIndex uiQueryBegin, nucSeqIndex uiQueryEnd)> fDo)
{//TODO: check git and continue here
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

void Bucketing::forEachNonBridgingPerfectMatch(std::shared_ptr<SegmentTreeInterval> pxNode, bool bAnchorOnly, std::shared_ptr<FM_Index> pxFM_index, std::shared_ptr<FM_Index> pxRev_FM_Index, std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence, std::shared_ptr<NucleotideSequence> pxQuerySeq, 
	std::function<void(std::shared_ptr<PerfectMatch>)> fDo)
{
	forEachNonBridgingHitOnTheRefSeq(
		pxNode, bAnchorOnly, pxFM_index, pxRev_FM_Index, pxRefSequence, pxQuerySeq,
		[&](nucSeqIndex ulIndexOnRefSeq, nucSeqIndex uiQuerryBegin, nucSeqIndex uiQuerryEnd)
		{
			fDo(std::shared_ptr<PerfectMatch>(new PerfectMatch(uiQuerryEnd - uiQuerryBegin, ulIndexOnRefSeq, uiQuerryBegin)));
		}//lambda
	);//for each
}//function

/* transfer the saved hits into the clustering
 * if DEBUG_CHECK_INTERVALS is activated the hits are verified before storing
*/
void Bucketing::saveHits(std::shared_ptr<SegmentTreeInterval> pxNode, std::shared_ptr<FM_Index> pxFM_index, std::shared_ptr<FM_Index> pxRev_FM_Index, std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence, std::shared_ptr<NucleotideSequence> pxQuerySeq, AnchorMatchList &rList)
{
	//bAnchorOnly = false since we also want to collet not maximally extended seeds
	forEachNonBridgingPerfectMatch(
		pxNode, false, pxFM_index, pxRev_FM_Index, pxRefSequence, pxQuerySeq,
		[&](std::shared_ptr<PerfectMatch> pxMatch)
		{
			rList.addMatch(std::shared_ptr<PerfectMatch>(pxMatch));
		}//lambda
	);//for each
}///function

/* transfer the anchors into the clustering
 * if DEBUG_CHECK_INTERVALS is activated the hits are verified before storing
*/
void Bucketing::saveAnchors(std::shared_ptr<SegmentTreeInterval> pxNode, std::shared_ptr<FM_Index> pxFM_index,  std::shared_ptr<FM_Index> pxRev_FM_Index, std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence, std::shared_ptr<NucleotideSequence> pxQuerySeq, AnchorMatchList &rList)
{
	//bAnchorOnly = true since we only want the maximally extended seeds
	forEachNonBridgingPerfectMatch(
		pxNode, true, pxFM_index, pxRev_FM_Index, pxRefSequence, pxQuerySeq,
		[&](std::shared_ptr<PerfectMatch> pxMatch)
		{
			rList.addAnchorSegment(std::shared_ptr<PerfectMatch>(pxMatch));
		}//lambda
	);//for each
}///function

std::shared_ptr<Container> Bucketing::execute(std::shared_ptr<Container> pInput)
{
	std::shared_ptr<ContainerVector> pCastedInput = std::static_pointer_cast<ContainerVector>(pInput);
	std::shared_ptr<SegmentTreeContainer> pSegments = std::static_pointer_cast<SegmentTreeContainer>(pCastedInput->vElements.at(0));
	std::shared_ptr<SegmentTreeContainer> pAnchors = std::static_pointer_cast<SegmentTreeContainer>(pCastedInput->vElements.at(1));
	std::shared_ptr<NucSeqContainer> pQuerrySeq = std::static_pointer_cast<NucSeqContainer>(pCastedInput->vElements.at(2));
	std::shared_ptr<PackContainer> pRefSeq = std::static_pointer_cast<PackContainer>(pCastedInput->vElements.at(3));
	std::shared_ptr<FM_IndexContainer> pFM_index = std::static_pointer_cast<FM_IndexContainer>(pCastedInput->vElements.at(4));
	std::shared_ptr<FM_IndexContainer> pFM_indexReversed = std::static_pointer_cast<FM_IndexContainer>(pCastedInput->vElements.at(5));

	AnchorMatchList xA(uiNumThreads, uiStripSize, pQuerrySeq->size(), pRefSeq->getUnpackedSize());

	/*
	*	extract all seeds from the segment tree intervals
	*	store them in buckets for easy pickup
	*/
	pSegments->pTree->forEach(
		[&](std::shared_ptr<SegmentTreeInterval> pxNode)
		{
			saveHits(pxNode, pFM_index->pIndex, pFM_indexReversed->pIndex, pRefSeq->pPack, pQuerrySeq->pSeq, xA);
		}//lambda
	);//forEach

	/*
	*	extract all anchor sequences
	*/
	pAnchors->pTree->forEach(
		[&](std::shared_ptr<SegmentTreeInterval> pxNode)
		{
			saveAnchors(pxNode, pFM_index->pIndex, pFM_indexReversed->pIndex, pRefSeq->pPack, pQuerrySeq->pSeq, xA);
		}//lambda
	);//forEach

	std::shared_ptr<StripOfConsiderationListContainer> pRet(new StripOfConsiderationListContainer());

	/*
	*	the actual work is hidden here:
	*		for each strip we pick up all hits lying in the respective buckets
	*/
	std::vector<std::shared_ptr<StripOfConsideration>> aStrips = xA.findAnchors();

	/*
	*	create the python container classes
	*/
	for(std::shared_ptr<StripOfConsideration> pStrip : aStrips)
	{
		pRet->push_back(std::shared_ptr<StripOfConsiderationContainer>(new StripOfConsiderationContainer(pStrip)));
	}//for

	return pRet;
}//function


void exportGraphicalMethod()
{
	//export the StripOfConsideration class
	boost::python::class_<
        StripOfConsiderationContainer, 
        boost::python::bases<Container>, 
        std::shared_ptr<StripOfConsiderationContainer>
    >("StripOfConsideration")
		.def("getScore", &StripOfConsiderationContainer::getValueOfContet)
		;

	//register a pointer to StripOfConsiderationContainer as return value to boost python
    boost::python::register_ptr_to_python< std::shared_ptr<StripOfConsiderationContainer> >();

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< std::shared_ptr<StripOfConsiderationContainer>, std::shared_ptr<Container> >(); 

	//export the StripOfConsiderationList class
	boost::python::class_<
        StripOfConsiderationListContainer, 
        boost::python::bases<Container>, 
        std::shared_ptr<StripOfConsiderationListContainer>
    >("StripOfConsiderationList")
		.def("size", &StripOfConsiderationListContainer::size)
		.def("at", &StripOfConsiderationListContainer::at)
		.def("append", &StripOfConsiderationListContainer::push_back)
		.def("remove", &StripOfConsiderationListContainer::remove)
		;
	//register a pointer to StripOfConsiderationListContainer as return value to boost python
    boost::python::register_ptr_to_python< std::shared_ptr<StripOfConsiderationListContainer> >();

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< std::shared_ptr<StripOfConsiderationListContainer>, std::shared_ptr<Container> >(); 

    //export the LineSweepContainer class
	//boost::python::class_<LineSweepContainer, boost::python::bases<Module>>("LineSweep")
	//	;
    //export the Bucketing class
	boost::python::class_<Bucketing, boost::python::bases<Module>>("Bucketing")
		;
}