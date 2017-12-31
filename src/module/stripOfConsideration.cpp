#include "module/stripOfConsideration.h"
using namespace libMABS;


ContainerVector StripOfConsideration::getInputType() const
{
    return ContainerVector
    {
        //all segments
        std::shared_ptr<Container>(new SegmentVector()),
        //the anchors
        std::shared_ptr<Container>(new Seeds()),
        //the query
        std::shared_ptr<Container>(new NucSeq()),
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
        std::shared_ptr<NucSeq> pxQuerySeq,
        std::function<void(Seed rxS)> fDo,
        nucSeqIndex addSize = 0
    )
{
    pVector->forEachSeed(
        pxFM_index, uiMaxAmbiguity, bSkipLongBWTIntervals,
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
        std::shared_ptr<ContainerVector> vpInput
    )
{
    std::shared_ptr<SegmentVector> pSegments = std::static_pointer_cast<SegmentVector>((*vpInput)[0]);
    std::shared_ptr<Seeds> pAnchors = std::static_pointer_cast<Seeds>((*vpInput)[1]);
    std::shared_ptr<NucSeq> pQuerySeq = 
        std::static_pointer_cast<NucSeq>((*vpInput)[2]);
    std::shared_ptr<Pack> pRefSeq = 
        std::static_pointer_cast<Pack>((*vpInput)[3]);
    std::shared_ptr<FMIndex> pFM_index = std::static_pointer_cast<FMIndex>((*vpInput)[4]);

    if (dMaxSeeds2 > 0 && pSegments->numSeeds(pFM_index, uiMaxAmbiguity) > pQuerySeq->length() * dMaxSeeds2)
        return std::shared_ptr<Seeds>(new Seeds());

    //extract the seeds
    std::vector<std::tuple<Seed, bool>> vSeeds;
    forEachNonBridgingSeed(
        pSegments, pFM_index, pRefSeq, pQuerySeq,
        [&](Seed xSeed)
        {
            vSeeds.push_back(std::make_tuple(xSeed, true));
        }//lambda
    );//for each

    //sort the seeds according to their initial positions
    std::sort(
        vSeeds.begin(), vSeeds.end(),
        [&]
        (const std::tuple<Seed, bool> a, const std::tuple<Seed, bool> b)
        {
            return getPositionForBucketing(pQuerySeq->length(), std::get<0>(a)) 
                    < getPositionForBucketing(pQuerySeq->length(), std::get<0>(b));
        }//lambda
    );//sort function call

    unsigned int uiAnchorIndex = 0;


    std::shared_ptr<ContainerVector> pRet(new ContainerVector(std::shared_ptr<Seeds>(new Seeds())));
    for(Seed& xAnchor : *pAnchors)
    {
        nucSeqIndex uiStart = getPositionForBucketing(pQuerySeq->length(), xAnchor) - uiStripSize/2;
        nucSeqIndex uiSize = uiStripSize;

        /*
        * FILTER START
        * we make sure that we can never have bridging strips
        */
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
        /*
        * FILTER END
        */

        std::shared_ptr<Seeds> pxNew(new Seeds());

        //binary search for the first element in range
        auto iterator = std::lower_bound(
            vSeeds.begin(), vSeeds.end(), uiStart,
            [&]
            (const std::tuple<Seed, bool> a, const nucSeqIndex uiStart)
            {
                return getPositionForBucketing(pQuerySeq->length(), std::get<0>(a)) <= uiStart;
            }//lambda
        );//binary search function call
        assert(getPositionForBucketing(pQuerySeq->length(), std::get<0>(*iterator)) >= uiStart);

        //save all seeds belonging into the strip of consideration
        while(
                iterator != vSeeds.end() &&
                getPositionForBucketing(pQuerySeq->length(), std::get<0>(*iterator)) < uiEnd
            )
        {
            /*
             * FILTER 2
             * we make sure that we don't collect the same seeds twice:
             */
            //if(std::get<1>(*iterator))
                pxNew->push_back(std::get<0>(*iterator));
            std::get<1>(*iterator) = false;
            ++iterator;
        }//while

        /*
         * FILTER 3
         */
        if(
                pxNew->size() < minSeeds && 
                pxNew->getScore() < minSeedLength * pQuerySeq->length() &&
                pSegments->numSeeds(pFM_index, uiMaxAmbiguity) > pQuerySeq->length() * dMaxSeeds
            )
            continue;


        /*
         * save some statistics about this strip
         */
        pxNew->xStats.index_of_strip = uiAnchorIndex++;
        pxNew->xStats.seed_coverage = pxNew->getScore();
        pxNew->xStats.num_seeds_in_strip = pxNew->size();
        pxNew->xStats.anchor_size = xAnchor.size();
        pxNew->xStats.anchor_ambiguity = xAnchor.uiAmbiguity;

        //save the strip of consideration
        pRet->push_back(pxNew);
    }//for
    return pRet;
}//function

void exportStripOfConsideration()
{
    //export the Bucketing class
    boost::python::class_<
            StripOfConsideration, 
            boost::python::bases<Module>, 
            std::shared_ptr<StripOfConsideration>
        >("StripOfConsideration")
        .def(boost::python::init<nucSeqIndex, unsigned int, unsigned int, float, float, float>())
        .def_readwrite("strip_size", &StripOfConsideration::uiStripSize)
        .def_readwrite("max_ambiguity", &StripOfConsideration::uiMaxAmbiguity)
        .def_readwrite("min_seeds", &StripOfConsideration::minSeeds)
        .def_readwrite("min_seed_length", &StripOfConsideration::minSeedLength)
        .def_readwrite("max_seeds", &StripOfConsideration::dMaxSeeds)
        .def_readwrite("max_seeds_2", &StripOfConsideration::dMaxSeeds2)
    ;

    boost::python::implicitly_convertible< 
        std::shared_ptr<StripOfConsideration>,
        std::shared_ptr<Module> 
    >();
}//function