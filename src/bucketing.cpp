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
	};
}//function

ContainerType Bucketing::getOutputType()
{
	return ContainerType::seedsVector;
}//function


void Bucketing::forEachNonBridgingSeed(
		std::shared_ptr<SegmentTreeInterval> pxNode,
		std::shared_ptr<FM_Index> pxFM_index,std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence,
		std::shared_ptr<NucleotideSequence> pxQuerySeq,
		std::function<void(Seed rxS)> fDo,
		nucSeqIndex addSize = 0
	)
{
	pxNode->forEachSeed(
		pxFM_index, uiMaxHitsPerInterval, bSkipLongBWTIntervals,
		[&](Seed xS)
		{
			// check if the match is bridging the forward/reverse strand 
			// or bridging between two chromosomes
			if ( pxRefSequence->bridingSubsection(
					//prevent negative index
					xS.start_ref() > addSize ? xS.start_ref() - addSize : 0,//from
					//prevent index larger than reference
					xS.end_ref() + addSize < pxFM_index->getRefSeqLength() ?
						xS.size() + addSize :
						pxFM_index->getRefSeqLength() - xS.start_ref() 
					) //to
				)
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
		std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence,
		std::shared_ptr<NucleotideSequence> pxQuerySeq,
		std::vector<SeedBucket>& raxSeedBuckets
	)
{
	forEachNonBridgingSeed(
		pxNode, pxFM_index, pxRefSequence, pxQuerySeq,
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
	std::shared_ptr<NucleotideSequence> pQuerySeq = 
		std::static_pointer_cast<NucleotideSequence>(vpInput[2]);
	std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pRefSeq = 
		std::static_pointer_cast<BWACompatiblePackedNucleotideSequencesCollection>(vpInput[3]);
	std::shared_ptr<FM_Index> pFM_index = std::static_pointer_cast<FM_Index>(vpInput[4]);

	
	std::vector<SeedBucket> axSeedBuckets((pQuerySeq->length() + pRefSeq->uiUnpackedSizeForwardPlusReverse())/uiStripSize + 1);


	/*
	*	extract all seeds from the segment tree intervals
	*	store them in buckets for easy pickup
	*/
	pSegments->forEach(
		[&](std::shared_ptr<SegmentTreeInterval> pxNode)
		{
			saveSeeds(pxNode, pFM_index, pRefSeq, pQuerySeq, axSeedBuckets);
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
				pxNode, pFM_index, pRefSeq, pQuerySeq,
				[&](Seed xAnchor)
				{
					nucSeqIndex uiStart = getPositionForBucketing(pQuerySeq->length(), xAnchor) - uiStripSize/2;
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
				},//lambda
				uiStripSize/2
			);//for each
		}//lambda
	);//forEach

	return pRet;
}//function

void exportBucketing()
{
    //export the Bucketing class
	boost::python::class_<Bucketing, boost::python::bases<CppModule>>(
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