#include "graphicalMethod.h"


std::vector<ContainerType> Bucketing::getInputType()
{
	return std::vector<ContainerType>
	{
		//all segments
		ContainerType::segmentList,
		//the anchors
		ContainerType::segmentList,
		//the querry
		ContainerType::nucSeq,
		//the reference
		ContainerType::packedNucSeq,
		//the forward fm_index
		ContainerType::fM_index,
		//the reversed fm_index
		ContainerType::fM_index
	};
}//function

std::vector<ContainerType> Bucketing::getOutputType()
{
	return std::vector<ContainerType>{ContainerType::stripOfConsiderationList};
}//function


void Bucketing::forEachNonBridgingSeed(
		std::shared_ptr<SegmentTreeInterval> pxNode,
		bool bAnchorOnly,
		std::shared_ptr<FM_Index> pxFM_index,
		std::shared_ptr<FM_Index> pxRev_FM_Index,std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence,
		std::shared_ptr<NucleotideSequence> pxQuerySeq,
		std::function<void(Seed rxS)> fDo
	)
{
	pxNode->forEachSeed(
		pxFM_index, pxRev_FM_Index, uiMaxHitsPerInterval, bSkipLongBWTIntervals, bAnchorOnly,
		[&](Seed xS)
		{
			//check if the match is bridging the forward/reverse strand or bridging between two chromosomes
			/* we have to make sure that the match does not start before or end after the reference sequence
			* this can happen since we can find parts on the end of the query at the very beginning of the reference or vis versa.
			* in this case we will replace the out of bounds index with 0 or the length of the reference sequence respectively.
			*/
			std::cout << "TODO: check bridging!" << std::endl;
			/*nucSeqIndex uiStart = xS.start_ref() - xS.start();
			if( xS.start_ref() < xS.start() )
				uiStart = 0;
			nucSeqIndex uiEnd = pxFM_index->getRefSeqLength() - ulIndexOnRefSeq;
			if(
					ulIndexOnRefSeq + pxQuerySeq->length() - xS.start()
						> 
					pxFM_index->getRefSeqLength()
				)
				uiEnd = pxFM_index->getRefSeqLength();
			if ( pxRefSequence->bridingSubsection(uiStart, uiEnd) )
			{
				//if so ignore this hit
				return;
			}//if*/
			fDo(xS);
		}//lambda
	);
}//function

/** transfer the saved hits into the clustering*/
void Bucketing::saveSeeds(
		std::shared_ptr<SegmentTreeInterval> pxNode,
		std::shared_ptr<FM_Index> pxFM_index,
		std::shared_ptr<FM_Index> pxRev_FM_Index,
		std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence,
		std::shared_ptr<NucleotideSequence> pxQuerySeq,
		std::vector<SeedBucket>& raxSeedBuckets
	)
{
	//bAnchorOnly = false since we also want to collet not maximally extended seeds
	forEachNonBridgingSeed(
		pxNode, false, pxFM_index, pxRev_FM_Index, pxRefSequence, pxQuerySeq,
		[&](Seed xSeed)
		{
			addSeed(pxQuerySeq->length(), xSeed, raxSeedBuckets);
		}//lambda
	);//for each
}///function

#if 0
/** transfer the anchors into the clustering
 * if DEBUG_CHECK_INTERVALS is activated the hits are verified before storing
*/
void Bucketing::saveAnchors(
		std::shared_ptr<SegmentTreeInterval> pxNode,
		std::shared_ptr<FM_Index> pxFM_index,  
		std::shared_ptr<FM_Index> pxRev_FM_Index, 
		std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence, 
		std::shared_ptr<NucleotideSequence> pxQuerySeq
	)
{
	//bAnchorOnly = true since we only want the maximally extended seeds
	forEachNonBridgingSeed(
		pxNode, true, pxFM_index, pxRev_FM_Index, pxRefSequence, pxQuerySeq,
		[&](const Seed xSeed)
		{
			addAnchor(xSeed);
		}//lambda
	);//for each
}///function
#endif

std::shared_ptr<Container> Bucketing::execute(
		std::vector<std::shared_ptr<Container>> vpInput
	)
{
	std::shared_ptr<SegmentTree> pSegments = std::static_pointer_cast<SegmentTree>(vpInput[0]);
	std::shared_ptr<SegmentTree> pAnchors = std::static_pointer_cast<SegmentTree>(vpInput[1]);
	std::shared_ptr<NucleotideSequence> pQuerrySeq = 
		std::static_pointer_cast<NucleotideSequence>(vpInput[2]);
	std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pRefSeq = 
		std::static_pointer_cast<BWACompatiblePackedNucleotideSequencesCollection>(vpInput[3]);
	std::shared_ptr<FM_Index> pFM_index = std::static_pointer_cast<FM_Index>(vpInput[4]);
	std::shared_ptr<FM_Index> pFM_indexReversed = std::static_pointer_cast<FM_Index>(vpInput[5]);

	
	std::vector<SeedBucket> axSeedBuckets((pQuerrySeq->length() + pRefSeq->uiUnpackedSizeForwardStrand)/uiStripSize);


	/*
	*	extract all seeds from the segment tree intervals
	*	store them in buckets for easy pickup
	*/
	pSegments->forEach(
		[&](std::shared_ptr<SegmentTreeInterval> pxNode)
		{
			saveSeeds(pxNode, pFM_index, pFM_indexReversed, pRefSeq, pQuerrySeq, axSeedBuckets);
		}//lambda
	);//forEach


	std::shared_ptr<StripOfConsiderationVector> pRet(new StripOfConsiderationVector());
	/*
	*	return one strip of consideration for each anchor
	*/
	pAnchors->forEach(
		[&](std::shared_ptr<SegmentTreeInterval> pxNode)
		{
			
			//bAnchorOnly = true since we want to collet maximally extended seeds
			forEachNonBridgingSeed(
				pxNode, true, pFM_index, pFM_indexReversed, pRefSeq, pQuerrySeq,
				[&](Seed xAnchor)
				{
					nucSeqIndex uiStart = getPositionForBucketing(pQuerrySeq->length(), xAnchor) - uiStripSize/2;
					std::shared_ptr<StripOfConsideration> pxNew(
							new StripOfConsideration(
								uiStart, 
								uiStripSize
							)
						);
					for (unsigned int uiC = pxNew->start() / uiStripSize - 1; uiC <= pxNew->end() / uiStripSize; uiC++)
					{
						pxNew->addElement(axSeedBuckets[uiC]);
					}//for
				}//lambda
			);//for each

			
		}//lambda
	);//forEach

	return pRet;
}//function
#if 0
std::vector<ContainerType> LineSweepContainer::getInputType()
{
	return std::vector<ContainerType>{
			//the querry
			ContainerType::nucSeq,
			//the reference
			ContainerType::packedNucSeq,
			//the stips of consideration
			ContainerType::stripOfConsiderationList,
		};
}//function

std::vector<ContainerType> LineSweepContainer::getOutputType()
{
	return std::vector<ContainerType>{ContainerType::stripOfConsideration};
}//function


std::shared_ptr<Container> LineSweepContainer::execute(
		std::vector<std::shared_ptr<Container>> vpInput
	)
{
	std::shared_ptr<NucleotideSequence> pQuerrySeq = 
		std::static_pointer_cast<NucleotideSequence>(vpInput[0]);
	std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pRefSeq =
		std::static_pointer_cast<BWACompatiblePackedNucleotideSequencesCollection>(vpInput[1]);
	std::shared_ptr<StripOfConsiderationVector> pStrips =
		std::static_pointer_cast<StripOfConsiderationVector>(vpInput[2]);

	GraphicalMethod xG(pRefSeq->uiUnpackedSizeForwardPlusReverse(), pQuerrySeq->length());

	/*
	*	extract the strips of consideration
	*/
	for(std::shared_ptr<StripOfConsideration> pContainer : pStrips->x)
	{
		xG.addStripOfConsideration(pContainer);
	}//for

	xG.smartProcess();

	return std::shared_ptr<StripOfConsiderationVector>(new StripOfConsiderationVector(xG.getNthBestBucket(0)));
}//function
#endif
void exportGraphicalMethod()
{
	//export the StripOfConsideration class
	boost::python::class_<
        StripOfConsideration, 
        boost::python::bases<Container>, 
        std::shared_ptr<StripOfConsideration>
    >(
			"StripOfConsideration",
			"Holds the matches close to a selected anchor match\n"
		);

	//register a pointer to StripOfConsideration as return value to boost python
    boost::python::register_ptr_to_python< std::shared_ptr<StripOfConsideration> >();

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< 
			std::shared_ptr<StripOfConsideration>, 
			std::shared_ptr<Container>
		>(); 

	//register return values of vectors of strips
	boost::python::class_<std::vector<std::shared_ptr<StripOfConsideration>>>("VecRetStrip")
		.def(boost::python::vector_indexing_suite<
				std::vector<std::shared_ptr<StripOfConsideration>>,
				/*
				*	true = noproxy this means that the content of the vector is already exposed by
				*	boost python. 
				*	if this is kept as false, StripOfConsideration would be exposed a second time.
				*	the two StripOfConsiderations would be different and not intercastable.
				*	=> keep this as true
				*/
				true
			>());

	//export the StripOfConsiderationVector class
	boost::python::class_<
			StripOfConsiderationVector, 
			boost::python::bases<Container>, 
			std::shared_ptr<StripOfConsiderationVector>
		>(
			"StripOfConsiderationVector",
			"	x: the vector holding the strips\n."
			"\n"
			"Contains multiple strips of consideration.\n"
		)
		.def_readwrite("x", &StripOfConsiderationVector::x);
	
	//register a pointer to StripOfConsideration as return value to boost python
	boost::python::register_ptr_to_python< std::shared_ptr<StripOfConsiderationVector> >();

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< 
			std::shared_ptr<StripOfConsiderationVector>, 
			std::shared_ptr<Container>
		>(); 
#if 0
    //export the LineSweepContainer class
	boost::python::class_<LineSweepContainer, boost::python::bases<Module>>(
			"LineSweep",
			"Uses linesweeping to remove contradicting "
			"matches within one strip of consideration.\n"
			"\n"
			"Execution:\n"
			"	Expects querry, ref, strip_vec as input.\n"
			"		querry: the querry as NucleotideSequence\n"
			"		ref: the reference seqeuence as Pack\n"
			"		strip_vec: the areas that shall be evaluated as StripOfConsiderationVector\n"
			"	returns strip_vec.\n"
			"		strip_vec: the evaluated areas\n"
		);
#endif
    //export the Bucketing class
	boost::python::class_<Bucketing, boost::python::bases<Module>>(
			"Bucketing",
			"Throws seeds into buckets in order to speed up the extraction "
			"of strips of consideration.\n"
			"\n"
			"	Execution:\n"
			"	Expects anchor, seeds, query, ref, ind, rev_ind as input.\n"
			"		anchor: the seeds that shall be used as anchors for the strips\n"
			"		seeds: all seeds\n"
			"		querry: the querry as NucleotideSequence\n"
			"		ref: the reference seqeuence as Pack\n"
			"		ind: the forward FM Index\n"
			"		rev_ind: the reversed FM Index\n"
			"	returns strip_vec.\n"
			"		strip_vec: a strip of consideration for each anchor "
			"as StripOfConsiderationVector\n"
		)
			.def_readwrite("strip_size", &Bucketing::uiStripSize)
			.def_readwrite("max_hits", &Bucketing::uiMaxHitsPerInterval)
			.def_readwrite("skip_long", &Bucketing::bSkipLongBWTIntervals);
}