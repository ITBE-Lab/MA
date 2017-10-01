#include "bucketing.h"



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
	return std::vector<ContainerType>{ContainerType::seedsVector};
}//function


void Bucketing::forEachNonBridgingSeed(
		std::shared_ptr<SegmentTreeInterval> pxNode,
		std::shared_ptr<FM_Index> pxFM_index,
		std::shared_ptr<FM_Index> pxRev_FM_Index,std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence,
		std::shared_ptr<NucleotideSequence> pxQuerySeq,
		std::function<void(Seed rxS)> fDo
	)
{
	pxNode->forEachSeed(
		pxFM_index, pxRev_FM_Index, uiMaxHitsPerInterval, bSkipLongBWTIntervals,
		[&](Seed xS)
		{
			//check if the match is bridging the forward/reverse strand or bridging between two chromosomes
			/* we have to make sure that the match does not start before or end after the reference sequence
			* this can happen since we can find parts on the end of the query at the very beginning of the reference or vis versa.
			* in this case we will replace the out of bounds index with 0 or the length of the reference sequence respectively.
			*/
			//TODO: check bridging!
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
	forEachNonBridgingSeed(
		pxNode, pxFM_index, pxRev_FM_Index, pxRefSequence, pxQuerySeq,
		[&](Seed xSeed)
		{
			addSeed(pxQuerySeq->length(), xSeed, raxSeedBuckets);
		}//lambda
	);//for each
}///function

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

	DEBUG(
		std::cout << "values of buckets: ";
		for(auto& bucket : axSeedBuckets)
		{
			std::cout << bucket.getValue() << " ";
		}//for
		std::cout << std::endl;
	)//DEBUG

	std::shared_ptr<SeedsVector> pRet(new SeedsVector());
	/*
	*	return one strip of consideration for each anchor
	*/
	pAnchors->forEach(
		[&](std::shared_ptr<SegmentTreeInterval> pxNode)
		{
			forEachNonBridgingSeed(
				pxNode, pFM_index, pFM_indexReversed, pRefSeq, pQuerrySeq,
				[&](Seed xAnchor)
				{
					nucSeqIndex uiStart = getPositionForBucketing(pQuerrySeq->length(), xAnchor) - uiStripSize/2;
					std::shared_ptr<Seeds> pxNew(new Seeds());
					for (
							unsigned int uiC = uiStart / uiStripSize - 1; 
							uiC <= (uiStart + uiStripSize) / uiStripSize; 
							uiC++
						)
					{
						for(const Seed& rS : axSeedBuckets[uiC].seeds())
							pxNew->push_back(rS);
					}//for
					pRet->push_back(pxNew);
				}//lambda
			);//for each
		}//lambda
	);//forEach

	return pRet;
}//function

void exportBucketing()
{
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
}//function