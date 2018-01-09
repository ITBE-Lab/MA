#include "module/stripOfConsideration.h"
using namespace libMABS;

extern int iGap;// = 16;
extern int iExtend;// = 1;
extern int iMatch;// = 8;
extern int iMissMatch;// = 2;


nucSeqIndex StripOfConsideration::getPositionForBucketing(nucSeqIndex uiQueryLength, const Seed xS)
{ 
    return xS.start_ref() + (uiQueryLength - xS.start()); 
}//function


nucSeqIndex StripOfConsideration::getStripSize(nucSeqIndex uiQueryLength)
{
    return (iMatch * uiQueryLength - iGap) / iExtend;
}//function

void StripOfConsideration::sort(std::vector<Seed>& vSeeds, nucSeqIndex qLen)
{
    //we need 34 bits max to express any index on any genome
    unsigned int max_bits_used = 34;
    unsigned int n = vSeeds.size();

    if( true /*2*max_bits_used*n / log2(n) > n*log2(n)*/ )
    {
        //quicksort is faster
        std::sort(
            vSeeds.begin(), vSeeds.end(),
            [&]
            (const Seed a, const Seed b)
            {
                return getPositionForBucketing(qLen, a) 
                        < getPositionForBucketing(qLen,b);
            }//lambda
        );//sort function call
    }//if
    else
    {
        //radix sort is faster than quicksort
        unsigned int amount_buckets = max_bits_used/log2(n);
        if(amount_buckets < 2)
            amount_buckets = 2;
        unsigned int iter = 0;
        std::vector<std::vector<Seed>> xBuckets1(
                amount_buckets, std::vector<Seed>()
            );
        std::vector<std::vector<Seed>> xBuckets2(
                amount_buckets, std::vector<Seed>()
            );
        for(auto& xSeed : vSeeds)
            xBuckets2[0].push_back(xSeed);
        auto pBucketsNow = &xBuckets1;
        auto pBucketsLast = &xBuckets2;

        //2^max_bits_used maximal possible genome size
        while(std::pow(amount_buckets,iter) <= std::pow(2,max_bits_used))
        {
            for(auto& xBucket : *pBucketsNow)
                xBucket.clear();
            for(auto& xBucket : *pBucketsLast)
                for(auto& xSeed : xBucket)
                {
                    nucSeqIndex index = getPositionForBucketing(qLen, xSeed);
                    index = (nucSeqIndex)( index / std::pow(amount_buckets,iter) ) % amount_buckets;
                    (*pBucketsNow)[index].push_back(xSeed);
                }//for
            iter++;
            //swap last and now
            auto pBucketsTemp = pBucketsNow;
            pBucketsNow = pBucketsLast;
            pBucketsLast = pBucketsTemp;
        }//while
        vSeeds.clear();
        DEBUG(
            nucSeqIndex last = 0;
        )//DEBUG
        for(auto& xBucket : *pBucketsLast)
            for(auto& xSeed : xBucket)
            {
                vSeeds.push_back(xSeed);
                DEBUG(
                    nucSeqIndex now = getPositionForBucketing(qLen, xSeed);
                    assert(now >= last);
                    last = now;
                )//DEBUG
            }//for
    }//else
}//function

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
                //returning true since we want to continue extracting seeds
                return true;
            }//if
            fDo(xS);
            //returning true since we want to continue extracting seeds
            return true;
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

    if (
            dMaxSeeds2 > 0 && 
            pSegments->numSeeds(pFM_index, uiMaxAmbiguity) > pQuerySeq->length() * dMaxSeeds2
        )
        return std::shared_ptr<Seeds>(new Seeds());

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
    sort(vSeeds, pQuerySeq->length());

    unsigned int uiAnchorIndex = 0;


    std::shared_ptr<ContainerVector> pRet(new ContainerVector(std::shared_ptr<Seeds>(new Seeds())));
    for(Seed& xAnchor : *pAnchors)
    {
        /*
         * This is the formula from the paper
         * computes the size required for the strip so that we collect all relevent seeds.
         */
        nucSeqIndex uiStripSize = getStripSize(pQuerySeq->length());


        nucSeqIndex uiStart = getPositionForBucketing(pQuerySeq->length(), xAnchor) - uiStripSize;
        nucSeqIndex uiSize = uiStripSize*2;

        /*
        * FILTER START
        * we make sure that we can never have bridging strips
        */
        if(uiStripSize > getPositionForBucketing(pQuerySeq->length(), xAnchor))
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
            (const Seed a, const nucSeqIndex uiStart)
            {
                return getPositionForBucketing(pQuerySeq->length(), a) <= uiStart;
            }//lambda
        );//binary search function call
        if(iterator != vSeeds.end())
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
        .def(boost::python::init<unsigned int, unsigned int, float, float, float>())
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