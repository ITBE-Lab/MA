#include "stripOfConsideration.h"



ContainerVector StripOfConsideration::getInputType() const
{
	return ContainerVector
	{
		//all segments
		std::shared_ptr<Container>(new SegmentVector()),
		//the anchors
		std::shared_ptr<Container>(new Seeds()),
		//the query
		std::shared_ptr<Container>(new NucleotideSequence()),
		//the reference
		std::shared_ptr<Container>(new Pack()),
		//the forward fm_index
		std::shared_ptr<Container>(new FMIndex()),
	};
}//function

std::shared_ptr<Container> StripOfConsideration::getOutputType() const
{
	return std::shared_ptr<Container>(new ContainerVector(
			std::shared_ptr<Container>(new Seeds())
		));
}//function


void StripOfConsideration::forEachNonBridgingSeed(
		std::shared_ptr<SegmentVector> pVector,
		std::shared_ptr<FMIndex> pxFM_index,std::shared_ptr<Pack> pxRefSequence,
		std::shared_ptr<NucleotideSequence> pxQuerySeq,
		std::function<void(Seed rxS)> fDo,
		nucSeqIndex addSize = 0
	)
{
	pVector->forEachSeed(
		pxFM_index, uiMaxHitsPerInterval, bSkipLongBWTIntervals,
		[&](Seed xS)
		{
			// check if the match is bridging the forward/reverse strand 
			// or bridging between two chromosomes
			if ( pxRefSequence->bridgingSubsection(
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


std::shared_ptr<Container> StripOfConsideration::execute(
		ContainerVector vpInput
	)
{
	std::shared_ptr<SegmentVector> pSegments = std::static_pointer_cast<SegmentVector>(vpInput[0]);
	std::shared_ptr<Seeds> pAnchors = std::static_pointer_cast<Seeds>(vpInput[1]);
	std::shared_ptr<NucleotideSequence> pQuerySeq = 
		std::static_pointer_cast<NucleotideSequence>(vpInput[2]);
	std::shared_ptr<Pack> pRefSeq = 
		std::static_pointer_cast<Pack>(vpInput[3]);
	std::shared_ptr<FMIndex> pFM_index = std::static_pointer_cast<FMIndex>(vpInput[4]);

	//extract the seeds
	std::vector<Seed> vSeeds;
	forEachNonBridgingSeed(
		pSegments, pFM_index, pRefSeq, pQuerySeq,
		[&](Seed xSeed)
		{
			vSeeds.push_back(xSeed);
		}//lambda
	);//for each

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
		if(pRefSeq->bridgingSubsection(uiStart, uiSize))
		{
			pRefSeq->unBridgeSubsection(uiStart, uiSize);
		}//if
		nucSeqIndex uiEnd = uiStart + uiSize;
		//2)
		//TODO: this has a squared complexity :(
		// -> might be saved by the fact that there are always very fey strips...?
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

		//save all seeds belonging into the strip of consideration
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

void exportStripOfConsideration()
{
    //export the Bucketing class
	boost::python::class_<
			StripOfConsideration, 
			boost::python::bases<CppModule>, 
            std::shared_ptr<StripOfConsideration>
		>(
        "StripOfConsideration"
    )
        .def_readwrite("strip_size", &StripOfConsideration::uiStripSize)
        .def_readwrite("max_hits", &StripOfConsideration::uiMaxHitsPerInterval)
        .def_readwrite("skip_long", &StripOfConsideration::bSkipLongBWTIntervals);

	boost::python::implicitly_convertible< 
		std::shared_ptr<StripOfConsideration>,
		std::shared_ptr<CppModule> 
	>();
}//function