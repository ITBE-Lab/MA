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

//@fixme compute the correct size and put in the formula here...
nucSeqIndex StripOfConsideration::getStripSize(nucSeqIndex uiQueryLength)
{
    return 1000;
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
            if ( !pxRefSequence->bridgingSubsection(
                    //prevent negative index
                    xS.start_ref() > addSize ? xS.start_ref() - addSize : 0,//from
                    //prevent index larger than reference
                    xS.end_ref() + addSize < pxFM_index->getRefSeqLength() ?
                        xS.size() + addSize :
                        pxFM_index->getRefSeqLength() - xS.start_ref() 
                    ) //to
                )
            {
                //if bridging ignore the hit
                fDo(xS);
            }//if
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
    std::shared_ptr<NucSeq> pQuerySeq = 
        std::static_pointer_cast<NucSeq>((*vpInput)[1]);
    std::shared_ptr<Pack> pRefSeq = 
        std::static_pointer_cast<Pack>((*vpInput)[2]);
    std::shared_ptr<FMIndex> pFM_index = std::static_pointer_cast<FMIndex>((*vpInput)[3]);

    nucSeqIndex uiQLen = pQuerySeq->length();
    
    /*
    * This is the formula from the paper
    * computes the size required for the strip so that we collect all relevent seeds.
    */
    nucSeqIndex uiStripSize = getStripSize(uiQLen);

    // FILTER
    if (
            dMaxSeeds > 0 &&
            pSegments->numSeeds(pFM_index, uiMaxAmbiguity) > uiQLen * dMaxSeeds
        )
        return std::shared_ptr<ContainerVector>(new ContainerVector(std::shared_ptr<Seeds>(new Seeds())));

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
    sort(vSeeds, uiQLen);

    //positions to remember the maxima
    std::vector<std::tuple<nucSeqIndex, std::vector<Seed>::iterator>> xMaxima;

    //find the SOC maxima
    nucSeqIndex uiCurrScore = 0;
    nucSeqIndex uiCurrEle = 0;
    std::vector<Seed>::iterator xStripStart = vSeeds.begin();
    std::vector<Seed>::iterator xStripEnd = vSeeds.begin();
    while(xStripStart != vSeeds.end() && xStripEnd != vSeeds.end())
    {
        //move xStripEnd forwards while it is closer to xStripStart than uiStripSize
        while(
            xStripEnd != vSeeds.end() &&
            getPositionForBucketing(uiQLen, *xStripStart) + uiStripSize 
            >= getPositionForBucketing(uiQLen, *xStripEnd))
        {
            //remember the additional score
            uiCurrScore += xStripEnd->getValue();
            //remember the additional element
            uiCurrEle++;
            //move the iterator forward
            xStripEnd++;
        }//while
        //here xStripEnd points one past the last element within the strip
        assert(uiCurrEle >= 1);
        //FILTER
        if(uiCurrEle > minSeeds || uiCurrScore > minSeedLength*uiQLen)
        {
            //check if we improved upon the last maxima while dealing with the same area
            if(
                !xMaxima.empty() && 
                getPositionForBucketing(uiQLen, *std::get<1>(xMaxima.back())) + uiStripSize
                >= getPositionForBucketing(uiQLen, *xStripStart))
            {
                if(std::get<0>(xMaxima.back()) < uiCurrScore)
                {
                    //if so we want to replace the old maxima
                    xMaxima.pop_back();
                    xMaxima.push_back(std::make_tuple(uiCurrScore, xStripStart));
                }//if
            }//if
            else
                //save the SOC
                xMaxima.push_back(std::make_tuple(uiCurrScore, xStripStart));
        }//if
        //move xStripStart one to the right (this will cause xStripEnd to be adjusted)
        assert(uiCurrScore >= xStripStart->getValue());
        uiCurrScore -= (xStripStart++)->getValue();
        uiCurrEle--;
    }//while

    // sort the SOC starting points according to the scores, 
    // so that we can extract the best SOC first
    std::sort(xMaxima.begin(), xMaxima.end(),
        []
        (std::tuple<nucSeqIndex, std::vector<Seed>::iterator> a, 
         std::tuple<nucSeqIndex, std::vector<Seed>::iterator> b)
        {
            return std::get<0>(a) > std::get<0>(b);
        }//lambda
    );//sort function call
    assert(xMaxima.size() <= 1 || std::get<0>(xMaxima.front()) >= std::get<0>(xMaxima.back()));

    //the collection of strips of consideration
    std::shared_ptr<ContainerVector> pRet(new ContainerVector(std::shared_ptr<Seeds>(new Seeds())));

    //extract the required amount of SOCs
    auto xCollect = xMaxima.begin();
    while(pRet->size() < numStrips && xCollect != xMaxima.end())
    {
        //the strip that shall be collected
        std::shared_ptr<Seeds> pSeeds(new Seeds());
        //iterator walking till the end of the strip that shall be collected
        auto xCollect2 = std::get<1>(*xCollect);
        nucSeqIndex end = getPositionForBucketing(uiQLen, *std::get<1>(*xCollect)) + uiStripSize;
        while(
            xCollect2 != vSeeds.end() &&
            end >= getPositionForBucketing(uiQLen, *xCollect2))
            //if the iterator is still within the strip add the seed and increment the iterator
            pSeeds->push_back(*(xCollect2++));
        //save the seed
        pRet->push_back(pSeeds);
        //move to the next strip
        xCollect++;
    }//while

    //return the strip collection
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
        .def(boost::python::init<unsigned int, unsigned int, unsigned int, float, float>())
        .def_readwrite("max_ambiguity", &StripOfConsideration::uiMaxAmbiguity)
        .def_readwrite("min_seeds", &StripOfConsideration::minSeeds)
        .def_readwrite("min_seed_length", &StripOfConsideration::minSeedLength)
        .def_readwrite("num_strips", &StripOfConsideration::numStrips)
        .def_readwrite("max_seeds", &StripOfConsideration::dMaxSeeds)
    ;

    boost::python::implicitly_convertible< 
        std::shared_ptr<StripOfConsideration>,
        std::shared_ptr<Module> 
    >();
}//function