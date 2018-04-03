/** 
 * @file stripOfConsideration.cpp
 * @brief Implements a Strip of Consideration.
 * @author Markus Schmidt
 */

#include "module/stripOfConsideration.h"
#include "container/minimizersHash.h"
using namespace libMA;

extern int iGap;
extern int iExtend;
extern int iMatch;
extern int iMissMatch;


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
        [&](const Seed& xS)
        {
            // check if the match is bridging the forward/reverse strand 
            // or bridging between two chromosomes
            if ( !pxRefSequence->bridgingSubsection(
                    //prevent negative index
                    xS.start_ref() > addSize ? xS.start_ref() - addSize : 0,//from
                    //prevent index larger than reference
                    xS.end_ref() + addSize < pxFM_index->getRefSeqLength() ?
                        xS.size() - 1 + addSize :
                        pxFM_index->getRefSeqLength() - xS.start_ref() - 1
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

/**
 * @brief used to determine more complex orders of SoCs 
 */
class SoCOrder
{
public:
    nucSeqIndex uiAccumulativeLength = 0;
    unsigned int uiSeedAmount = 0;

    inline void operator+=(const Seed& rS)
    {
        uiSeedAmount++;
        uiAccumulativeLength += rS.getValue();
    }//operator

    inline void operator-=(const Seed& rS)
    {
        assert(uiSeedAmount > 0);
        uiSeedAmount--;
        assert(uiAccumulativeLength >= rS.getValue());
        uiAccumulativeLength -= rS.getValue();
    }//operator

    inline bool operator<(const SoCOrder& rOther) const
    {
        if(uiAccumulativeLength == rOther.uiAccumulativeLength)
            return uiSeedAmount > rOther.uiSeedAmount;
        return uiAccumulativeLength < rOther.uiAccumulativeLength;
    }//operator

    inline void operator=(const SoCOrder& rOther)
    {
        uiAccumulativeLength = rOther.uiAccumulativeLength;
        uiSeedAmount = rOther.uiSeedAmount;
    }//operator
}; //class

std::shared_ptr<Container> StripOfConsideration::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
    const std::shared_ptr<SegmentVector>& pSegments = std::static_pointer_cast<SegmentVector>((*vpInput)[0]);
    const std::shared_ptr<NucSeq>& pQuerySeq = 
        std::static_pointer_cast<NucSeq>((*vpInput)[1]);
    const std::shared_ptr<Pack>& pRefSeq = 
        std::static_pointer_cast<Pack>((*vpInput)[2]);
    const std::shared_ptr<FMIndex>& pFM_index = std::static_pointer_cast<FMIndex>((*vpInput)[3]);

    const nucSeqIndex uiQLen = pQuerySeq->length();
    
    /*
    * This is the formula from the paper
    * computes the size required for the strip so that we collect all relevent seeds.
    */
    const nucSeqIndex uiStripSize = this->getStripSize(uiQLen, iMatch, iExtend, iGap);

    //extract the seeds
    std::vector<Seed> vSeeds;
    vSeeds.reserve(pSegments->numSeeds(uiMaxAmbiguity));
    forEachNonBridgingSeed(
        pSegments, pFM_index, pRefSeq, pQuerySeq,
        [&](const Seed& xSeed)
        {
            vSeeds.push_back(xSeed);
        }//lambda
    );//for each

    //make sure that we return at least an empty seed set if nothing else
    if(vSeeds.empty())
        return std::shared_ptr<ContainerVector>(
            new ContainerVector {std::shared_ptr<Seeds>(new Seeds())}
        );

    //sort the seeds according to their initial positions
    std::sort(
        vSeeds.begin(), vSeeds.end(),
        [&]
        (const Seed &a, const Seed &b)
        {
            return getPositionForBucketing(uiQLen, a) 
                    < getPositionForBucketing(uiQLen,b);
        }//lambda
    );//sort function call

    //positions to remember the maxima
    std::vector<std::pair<SoCOrder, std::vector<Seed>::iterator>> vMaxima;

    //find the SOC maxima
    SoCOrder xCurrScore;
    std::vector<Seed>::iterator xStripStart = vSeeds.begin();
    std::vector<Seed>::iterator xStripEnd = vSeeds.begin();
    while(xStripEnd != vSeeds.end())
    {
        //move xStripEnd forwards while it is closer to xStripStart than uiStripSize
        nucSeqIndex uiCurrSize = 0;
        while(
            xStripEnd != vSeeds.end() &&
            getPositionForBucketing(uiQLen, *xStripStart) + uiStripSize 
            >= getPositionForBucketing(uiQLen, *xStripEnd))
        {
            //remember the additional score
            xCurrScore += *xStripEnd;
            // compute the current SOC size
            uiCurrSize = xStripStart->start_ref() - xStripEnd->start_ref();
            // carefull here we might have seeds in the wrong order since we sorted by r - q not r.
            if(xStripEnd->start_ref() > xStripStart->start_ref())
                uiCurrSize = xStripEnd->start_ref() - xStripStart->start_ref();
            //move the iterator forward
            xStripEnd++;
        }//while
        //here xStripEnd points one past the last element within the strip
        int64_t iDummy;
        //FILTER
        if(!pRefSeq->bridgingSubsection(xStripStart->start_ref(), uiCurrSize, iDummy))
        {
            //check if we improved upon the last maxima while dealing with the same area
            if(
                !vMaxima.empty() && 
                getPositionForBucketing(uiQLen, *vMaxima.back().second) + uiStripSize
                >= getPositionForBucketing(uiQLen, *xStripStart))
            {
                if(vMaxima.back().first < xCurrScore)
                {
                    //if so we want to replace the old maxima
                    vMaxima.pop_back();
                    vMaxima.push_back(std::make_pair(xCurrScore, xStripStart));
                }//if
            }//if
            else
                //save the SOC
                vMaxima.push_back(std::make_pair(xCurrScore, xStripStart));
        }//if
        //move xStripStart one to the right (this will cause xStripEnd to be adjusted)
        xCurrScore -= *(xStripStart++);
    }//while


    // make max heap from the SOC starting points according to the scores, 
    // so that we can extract the best SOC first
    auto vHeapOrder = 
        []
        (
            const std::pair<SoCOrder, std::vector<Seed>::iterator>& rA,
            const std::pair<SoCOrder, std::vector<Seed>::iterator>& rB
        )
        {
            return rA.first < rB.first;
        }//lambda
    ;
    std::make_heap(vMaxima.begin(), vMaxima.end(), vHeapOrder);//make heap function call

    //the collection of strips of consideration
    std::shared_ptr<ContainerVector> pRet(new ContainerVector(std::shared_ptr<Seeds>(new Seeds())));

    /*
     * @todo: move this into the SoC loop!!! (should save an insane amount of runtime)
     * if the best SoC quality is lower than fGiveUp * uiQLen we give up the entire process here
     * fGiveUp = 0 disables this.
     */
    if(fGiveUp != 0 && vMaxima.front().first.uiAccumulativeLength < fGiveUp * uiQLen)
        return pRet;

    unsigned int uiSoCIndex = 0;
    //extract the required amount of SOCs
    while(pRet->size() < numStrips && !vMaxima.empty())
    {
        //the strip that shall be collected
        std::shared_ptr<Seeds> pSeeds(new Seeds());
        // get the expected amount of seeds in the SoC from the order class and reserve memory
        pSeeds->reserve(vMaxima.front().first.uiSeedAmount);
        //iterator walking till the end of the strip that shall be collected
        auto xCollect2 = vMaxima.front().second;
        //save SoC index
        pSeeds->xStats.index_of_strip = uiSoCIndex++;
        // all these things are not used at the moment...
        pSeeds->xStats.uiInitialQueryBegin = xCollect2->start();
        pSeeds->xStats.uiInitialRefBegin = xCollect2->start_ref();
        pSeeds->xStats.uiInitialQueryEnd = xCollect2->end();
        pSeeds->xStats.uiInitialRefEnd = xCollect2->end_ref();
        nucSeqIndex end = getPositionForBucketing(uiQLen, *xCollect2) + uiStripSize;
        while(
            xCollect2 != vSeeds.end() &&
            end >= getPositionForBucketing(uiQLen, *xCollect2))
        {
            // save the beginning and end of the SoC
            // all these things are not used at the moment...
            if(xCollect2->start() < pSeeds->xStats.uiInitialQueryBegin)
                pSeeds->xStats.uiInitialQueryBegin = xCollect2->start();
            if(xCollect2->start_ref() < pSeeds->xStats.uiInitialRefBegin)
                pSeeds->xStats.uiInitialRefBegin = xCollect2->start_ref();
            if(xCollect2->end() > pSeeds->xStats.uiInitialQueryEnd)
                pSeeds->xStats.uiInitialQueryEnd = xCollect2->end();
            if(xCollect2->end_ref() > pSeeds->xStats.uiInitialRefEnd)
                pSeeds->xStats.uiInitialRefEnd = xCollect2->end_ref();
            pSeeds->xStats.num_seeds_in_strip++;
            assert(xCollect2->start() <= xCollect2->end());
            assert(xCollect2->end() <= pQuerySeq->length());
            //if the iterator is still within the strip add the seed and increment the iterator
            pSeeds->push_back(*(xCollect2++));
        }//while
        //save the seed
        pRet->push_back(pSeeds);
        //move to the next strip
        std::pop_heap (vMaxima.begin(), vMaxima.end(), vHeapOrder); vMaxima.pop_back();
    }//while

    //make sure that we return at least an empty seed set if nothing else
    if(pRet->size() == 0)
        pRet->push_back(std::shared_ptr<Seeds>(new Seeds()));

    //return the strip collection
    return pRet;
}//function

nucSeqIndex StripOfConsideration2::getPositionForBucketing(nucSeqIndex uiQueryLength, const Seed xS)
{ 
    return xS.start_ref() + (uiQueryLength - xS.start()); 
}//function

nucSeqIndex StripOfConsideration2::getStripSize(nucSeqIndex uiQueryLength)
{
    return (iMatch * uiQueryLength - iGap) / iExtend;
}//function

ContainerVector StripOfConsideration2::getInputType() const
{
    return ContainerVector
    {
        //the seeds
        std::shared_ptr<Container>(new Seeds()),
        //the query
        std::shared_ptr<Container>(new NucSeq()),
        //the reference
        std::shared_ptr<Container>(new Pack()),
    };
}//function

std::shared_ptr<Container> StripOfConsideration2::getOutputType() const
{
    return std::shared_ptr<Container>(new ContainerVector(
            std::shared_ptr<Container>(new Seeds())
        ));
}//function

std::shared_ptr<Container> StripOfConsideration2::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
    const std::shared_ptr<Seeds>& pSeeds = std::static_pointer_cast<Seeds>((*vpInput)[0]);
    const std::shared_ptr<NucSeq>& pQuerySeq = 
        std::static_pointer_cast<NucSeq>((*vpInput)[1]);
    const std::shared_ptr<Pack>& pRefSeq = 
        std::static_pointer_cast<Pack>((*vpInput)[2]);

    nucSeqIndex uiQLen = pQuerySeq->length();
    
    /*
    * This is the formula from the paper
    * computes the size required for the strip so that we collect all relevent seeds.
    */
    const nucSeqIndex uiStripSize = getStripSize(uiQLen);

    DEBUG(
        //verify that all seeds are correct
        for(auto xSeed : *pSeeds)
        {
            auto qSeg = pQuerySeq->fromTo(xSeed.start(), xSeed.end());
            auto rSeg = pRefSeq->vExtract(xSeed.start_ref(), xSeed.end_ref());
            if(qSeg != rSeg->toString())
            {
                std::cout
                    << "faulty seed: "
                    << qSeg
                    << " != "
                    << rSeg->toString() 
                    << " on rev comp strand: "
                    << (pRefSeq->bPositionIsOnReversStrand(xSeed.start_ref()) ? "true " : "false")
                    << " indices match: "
                    //watch out this is statically set to 15 in the moment
                    << ( Minimizer<15>(pQuerySeq,xSeed.start()) == Minimizer<15>(rSeg,0) ? "true " : "false")
                    << std::endl;
                //assert(false);
            }//if
        }//for
    )//DEBUG


    //make sure that we return at least an empty seed set if nothing else
    if(pSeeds->empty())
        return std::shared_ptr<ContainerVector>(
            new ContainerVector {std::shared_ptr<Seeds>(new Seeds())}
        );

    //positions to remember the maxima
    std::vector<std::tuple<nucSeqIndex, std::vector<Seed>::iterator, unsigned long>> vMaxima;

    //find the SOC maxima
    nucSeqIndex uiCurrScore = 0;
    unsigned long uiCurrEle = 0;
    std::vector<Seed>::iterator xStripStart = pSeeds->begin();
    std::vector<Seed>::iterator xStripEnd = pSeeds->begin();
    while(xStripEnd != pSeeds->end())
    {
        //move xStripEnd forwards while it is closer to xStripStart than uiStripSize
        nucSeqIndex uiCurrSize = 0;
        while(
            xStripEnd != pSeeds->end() &&
            getPositionForBucketing(uiQLen, *xStripStart) + uiStripSize 
            >= getPositionForBucketing(uiQLen, *xStripEnd))
        {
            //remember the additional score
            uiCurrScore += xStripEnd->getValue();
            // compute the current SOC size
            uiCurrSize = xStripStart->start_ref() - xStripEnd->start_ref();
            // carefull here we might have seeds in the wrong order since we sorted by r - q not r.
            if(xStripEnd->start_ref() > xStripStart->start_ref())
                uiCurrSize = xStripEnd->start_ref() - xStripStart->start_ref();
            //remember the additional element
            uiCurrEle++;
            //move the iterator forward
            xStripEnd++;
        }//while
        //here xStripEnd points one past the last element within the strip
        assert(uiCurrEle >= 1);
        int64_t iDummy;
        //FILTER
        if(!pRefSeq->bridgingSubsection(xStripStart->start_ref(), uiCurrSize, iDummy))
        {
            //check if we improved upon the last maxima while dealing with the same area
            if(
                !vMaxima.empty() && 
                getPositionForBucketing(uiQLen, *std::get<1>(vMaxima.back())) + uiStripSize
                >= getPositionForBucketing(uiQLen, *xStripStart))
            {
                if(std::get<0>(vMaxima.back()) < uiCurrScore)
                {
                    //if so we want to replace the old maxima
                    vMaxima.pop_back();
                    vMaxima.push_back(std::make_tuple(uiCurrScore, xStripStart, uiCurrEle));
                }//if
            }//if
            else
                //save the SOC
                vMaxima.push_back(std::make_tuple(uiCurrScore, xStripStart, uiCurrEle));
        }//if
        //move xStripStart one to the right (this will cause xStripEnd to be adjusted)
        assert(uiCurrScore >= xStripStart->getValue());
        uiCurrScore -= (xStripStart++)->getValue();
        uiCurrEle--;
    }//while

    // sort the SOC starting points according to the scores, 
    // so that we can extract the best SOC first
    std::sort(vMaxima.begin(), vMaxima.end(),
        []
        (std::tuple<nucSeqIndex, std::vector<Seed>::iterator, unsigned long> a, 
         std::tuple<nucSeqIndex, std::vector<Seed>::iterator, unsigned long> b)
        {
            if(std::get<0>(a) == std::get<0>(b))
                return std::get<2>(a) < std::get<2>(b);
            return std::get<0>(a) > std::get<0>(b);
        }//lambda
    );//sort function call
    assert(vMaxima.size() <= 1 || std::get<0>(vMaxima.front()) >= std::get<0>(vMaxima.back()));

    //the collection of strips of consideration
    std::shared_ptr<ContainerVector> pRet(new ContainerVector(std::shared_ptr<Seeds>(new Seeds())));

    //extract the required amount of SOCs
    auto xCollect = vMaxima.begin();
    while(pRet->size() < numStrips && xCollect != vMaxima.end())
    {
        //the strip that shall be collected
        std::shared_ptr<Seeds> pSeedsNew(new Seeds());
        //iterator walking till the end of the strip that shall be collected
        auto xCollect2 = std::get<1>(*xCollect);
        nucSeqIndex end = getPositionForBucketing(uiQLen, *std::get<1>(*xCollect)) + uiStripSize;
        while(
            xCollect2 != pSeeds->end() &&
            end >= getPositionForBucketing(uiQLen, *xCollect2))
        {
            assert(xCollect2->start() <= xCollect2->end());
            if(xCollect2->end() > pQuerySeq->length())
                std::cout << xCollect2->end() << std::endl;
            assert(xCollect2->end() <= pQuerySeq->length());
            //if the iterator is still within the strip add the seed and increment the iterator
            pSeedsNew->push_back(*(xCollect2++));
        }//while
        //save the seed
        pRet->push_back(pSeedsNew);
        //move to the next strip
        xCollect++;
    }//while

    //make sure that we return at least an empty seed set if nothing else
    if(pRet->size() == 0)
        pRet->push_back(std::shared_ptr<Seeds>(new Seeds()));

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
        .def(boost::python::init<unsigned int, unsigned int, float, float>())
        .def_readwrite("max_ambiguity", &StripOfConsideration::uiMaxAmbiguity)
        .def_readwrite("num_strips", &StripOfConsideration::numStrips)
        .def_readonly("min_score", &StripOfConsideration::fScoreMinimum)
    ;

    boost::python::implicitly_convertible< 
        std::shared_ptr<StripOfConsideration>,
        std::shared_ptr<Module> 
    >();
    //export the Bucketing class
    boost::python::class_<
            StripOfConsideration2, 
            boost::python::bases<Module>, 
            std::shared_ptr<StripOfConsideration2>
        >("StripOfConsideration2")
        .def(boost::python::init<unsigned int, unsigned int>())
        .def_readwrite("max_ambiguity", &StripOfConsideration2::uiMaxAmbiguity)
        .def_readwrite("num_strips", &StripOfConsideration2::numStrips)
    ;

    boost::python::implicitly_convertible< 
        std::shared_ptr<StripOfConsideration2>,
        std::shared_ptr<Module> 
    >();
}//function