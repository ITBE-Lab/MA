#include "bucketing.h"



std::vector<std::shared_ptr<Container>> Bucketing::getInputType() const
{
	return std::vector<std::shared_ptr<Container>>
	{
		//all segments
		std::shared_ptr<Container>(new SegmentTree()),
		//the anchors
		std::shared_ptr<Container>(new Seeds()),
		//the query
		std::shared_ptr<Container>(new NucleotideSequence()),
		//the reference
		std::shared_ptr<Container>(new BWACompatiblePackedNucleotideSequencesCollection()),
		//the forward fm_index
		std::shared_ptr<Container>(new FM_Index()),
	};
}//function

std::shared_ptr<Container> Bucketing::getOutputType() const
{
	return std::shared_ptr<Container>(new ContainerVector(
			std::shared_ptr<Container>(new Seeds())
		));
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
	std::shared_ptr<Seeds> pAnchors = std::static_pointer_cast<Seeds>(vpInput[1]);
	std::shared_ptr<NucleotideSequence> pQuerySeq = 
		std::static_pointer_cast<NucleotideSequence>(vpInput[2]);
	std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pRefSeq = 
		std::static_pointer_cast<BWACompatiblePackedNucleotideSequencesCollection>(vpInput[3]);
	std::shared_ptr<FM_Index> pFM_index = std::static_pointer_cast<FM_Index>(vpInput[4]);

	std::vector<Seed> vSeeds;
	pSegments->forEach(
		[&](std::shared_ptr<SegmentTreeInterval> pxNode)
		{
			forEachNonBridgingSeed(
				pxNode, pFM_index, pRefSeq, pQuerySeq,
				[&](Seed xSeed)
				{
					vSeeds.push_back(xSeed);
				}//lambda
			);//for each
		}//lambda
	);//forEach

	//sort the seeds according to their initial positions
	std::sort(
		vSeeds.begin(), vSeeds.end(),
		[&]
		(const Seed a, const Seed b)
		{
			return getPositionForBucketing(pQuerySeq->length(), a) 
					< getPositionForBucketing(pQuerySeq->length(), b);
		}//lambda
	);//sort function call

	//used to make sure we don't collect the same area twice
	std::vector<std::tuple<nucSeqIndex, nucSeqIndex>> collectedIntervals;

	std::shared_ptr<ContainerVector> pRet(new ContainerVector(std::shared_ptr<Seeds>(new Seeds())));
	for(Seed& xAnchor : *pAnchors)
	{
		nucSeqIndex uiStart = getPositionForBucketing(pQuerySeq->length(), xAnchor) - uiStripSize/2;
		nucSeqIndex uiSize = uiStripSize;

		/*
			* FILTER START
			*
			* 1)	we make sure that we can never have bridging strips
			*
			* 2) 	we filter out anchors that are to close to each other
			*		since we do not want to work on the same area twice
			*/
		//1)
		if(uiStripSize/2 > getPositionForBucketing(pQuerySeq->length(), xAnchor))
			uiStart = 0;
		if(uiStart >= pRefSeq->uiUnpackedSizeForwardPlusReverse())
			uiStart = pRefSeq->uiUnpackedSizeForwardPlusReverse() - 1;
		if(uiStart + uiSize >= pRefSeq->uiUnpackedSizeForwardPlusReverse())
			uiSize = pRefSeq->uiUnpackedSizeForwardPlusReverse() - uiStart - 1;
		if(pRefSeq->bridingSubsection(uiStart, uiSize))
		{
			pRefSeq->unBridgeSubsection(uiStart, uiSize);
		}//if
		nucSeqIndex uiEnd = uiStart + uiSize;
		//2)
		for(std::tuple<nucSeqIndex, nucSeqIndex> intv : collectedIntervals)
		{
			if(std::get<0>(intv) >= uiEnd)
				continue;
			if(std::get<1>(intv) <= uiStart)
				continue;
		}//for
		collectedIntervals.push_back(std::make_tuple(uiStart, uiEnd));
		/*
			* FILTER END
			*/

		std::shared_ptr<Seeds> pxNew(new Seeds());

		//binary search for the first element in range
		auto iterator = std::lower_bound(
			vSeeds.begin(), vSeeds.end(), uiStart,
			[&]
			(const Seed a, const nucSeqIndex uiStart)
			{
				return getPositionForBucketing(pQuerySeq->length(), a) < uiStart;
			}//lambda
		);//binary search function call
		assert(getPositionForBucketing(pQuerySeq->length(), *iterator) >= uiStart);

		while(
				iterator != vSeeds.end() &&
				getPositionForBucketing(pQuerySeq->length(), *iterator) < uiEnd
			)
		{
			pxNew->push_back(*iterator);
			++iterator;
		}//while

		pRet->push_back(pxNew);
	}//for
	return pRet;
}//function

void exportBucketing()
{
    //export the Bucketing class
	boost::python::class_<
			Bucketing, 
			boost::python::bases<CppModule>, 
            std::shared_ptr<Bucketing>
		>(
        "Bucketing",
        "Throws seeds into buckets in order to speed up the extraction "
        "of strips of consideration.\n"
        "\n"
        "	Execution:\n"
        "	Expects anchor, seeds, query, ref, ind, rev_ind as input.\n"
        "		anchor: the seeds that shall be used as anchors for the strips\n"
        "		seeds: all seeds\n"
        "		query: the query as NucleotideSequence\n"
        "		ref: the reference sequence as Pack\n"
        "		ind: the forward FM Index\n"
        "		rev_ind: the reversed FM Index\n"
        "	returns strip_vec.\n"
        "		strip_vec: a strip of consideration for each anchor "
        "as StripOfConsiderationVector\n"
    )
        .def_readwrite("strip_size", &Bucketing::uiStripSize)
        .def_readwrite("max_hits", &Bucketing::uiMaxHitsPerInterval)
        .def_readwrite("skip_long", &Bucketing::bSkipLongBWTIntervals);

	boost::python::implicitly_convertible< 
		std::shared_ptr<Bucketing>,
		std::shared_ptr<CppModule> 
	>();
}//function