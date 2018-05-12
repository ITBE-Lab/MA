/** 
 * @file stripOfConsideration.cpp
 * @brief Implements a Strip of Consideration.
 * @author Markus Schmidt
 */

#include "module/stripOfConsideration.h"
#include "util/system.h"
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
    return std::shared_ptr<Container>(new SoCPriorityQueue());
}//function

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
    auto pSeeds = std::make_shared<Seeds>();
    pSeeds->xStats.sName = pQuerySeq->sName;
    // rough estimate of how many seeds we will have 
    // (trying to avoid multiple allocations)
    pSeeds->reserve(pSegments->size() * 3);
    
#if MEASURE_DURATIONS == ( 1 ) 
    pSegments->fExtraction += metaMeasureDuration ( [&] () 
    {
#endif
            emplaceAllNonBridgingSeed(
                    *pSegments, // Segment vector (outcome of seeding)
                    *pFM_index,
                    *pRefSeq,
                    *pSeeds,
                    uiQLen
                );//emplace all function call
#if MEASURE_DURATIONS == ( 1 )
        }// lambda
    ).count() * 1000 ;// metaMeasureDuration function call
#endif
    

    //make sure that we return at least an SoC set
    if(pSeeds->empty())
        return std::make_shared<SoCPriorityQueue>(uiStripSize, pSeeds);

#if MEASURE_DURATIONS == ( 1 )
    pSegments->fSorting += metaMeasureDuration(
        [&]
        ()
        {
#endif


            //sort the seeds according to their initial positions
            std::sort(
                pSeeds->begin(), pSeeds->end(),
                [&]
                (const Seed &a, const Seed &b)
                {
        #if DELTA_CACHE == ( 1 )
                    return a.uiDelta < b.uiDelta;
        #else
                    return    getPositionForBucketing( uiQLen, a ) 
                            < getPositionForBucketing( uiQLen, b );
        #endif
                }//lambda
            );//sort function call


#if MEASURE_DURATIONS == ( 1 )
        }// lambda
    ).count() * 1000 ;// metaMeasureDuration function call
#endif

    //positions to remember the maxima
    auto pSoCs = std::make_shared<SoCPriorityQueue>(uiStripSize, pSeeds);

#if MEASURE_DURATIONS == ( 1 )
    pSegments->fLinesweep += metaMeasureDuration(
        [&]
        ()
        {
#endif


            //find the SOC maxima
            SoCOrder xCurrScore;
            std::vector<Seed>::iterator xStripStart = pSeeds->begin();
            std::vector<Seed>::iterator xStripEnd = pSeeds->begin();
            while(xStripEnd != pSeeds->end() && xStripStart != pSeeds->end())
            {
                //move xStripEnd forwards while it is closer to xStripStart than uiStripSize
                nucSeqIndex uiCurrSize = 0;
                while(
                    xStripEnd != pSeeds->end() &&
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


                DEBUG(
                    pSoCs->vScores.push_back(std::make_pair(
                        xCurrScore.uiAccumulativeLength, 
                        xStripStart->start_ref())
                    );
                ) // DEBUG


                //here xStripEnd points one past the last element within the strip
                int64_t iDummy;

                //FILTER
                /*
                * if the SoC quality is lower than fGiveUp * uiQLen we do not consider this SoC at all
                * fGiveUp = 0 disables this.
                */
                if(
                        ( fGiveUp == 0 || xCurrScore.uiAccumulativeLength >= fGiveUp * uiQLen ) && 
                        !pRefSeq->bridgingSubsection(xStripStart->start_ref(), uiCurrSize, iDummy)
                    )
                    pSoCs->push_back_no_overlap(
                        xCurrScore,
                        xStripStart,
                        xStripEnd,
                        getPositionForBucketing(uiQLen, *xStripStart),
                        getPositionForBucketing( uiQLen, *(xStripEnd-1) )
                    );
                // move xStripStart one to the right (this will cause xStripEnd to be adjusted)
                xCurrScore -= *(xStripStart++);
            }// while


            // make a max heap from the SOC starting points according to the scores, 
            // so that we can extract the best SOC first
            pSoCs->make_heap();


#if MEASURE_DURATIONS == ( 1 )
        }// lambda
    ).count() * 1000 ;// metaMeasureDuration function call
#endif

    //return the strip collection
    return pSoCs;

}//function

#ifdef WITH_PYTHON
void exportStripOfConsideration()
{
    //export the Bucketing class
    boost::python::class_<
            StripOfConsideration, 
            boost::python::bases<Module>, 
            std::shared_ptr<StripOfConsideration>
        >("StripOfConsideration")
        .def(boost::python::init<float>())
        .def_readwrite("max_ambiguity", &StripOfConsideration::uiMaxAmbiguity)
        .def_readwrite("min_score", &StripOfConsideration::fScoreMinimum)
        .def_readwrite("skip_long_bwt_intervals", &StripOfConsideration::bSkipLongBWTIntervals)
    ;

    boost::python::implicitly_convertible< 
        std::shared_ptr<StripOfConsideration>,
        std::shared_ptr<Module> 
    >();

}//function
#endif
