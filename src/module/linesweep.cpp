/** 
 * @file linesweep.cpp
 * @author Markus Schmidt
 */
#include "module/linesweep.h"
#if USE_RANSAC == 1
    #include "sample_consensus/test_ransac.h"
#endif
using namespace libMA;

extern int iGap;
extern int iExtend;
extern int iMatch;
extern int iMissMatch;
extern nucSeqIndex uiMaxGapArea;

ContainerVector LinearLineSweep::getInputType() const
{
    return ContainerVector{
            std::make_shared<SoCPriorityQueue>(),
            std::make_shared<NucSeq>()
        };
}// function

std::shared_ptr<Container> LinearLineSweep::getOutputType() const
{
    return std::make_shared<ContainerVector>(std::make_shared<Seeds>());
}// function

std::shared_ptr<std::vector<std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>>>
    LinearLineSweep::linesweep(
        std::shared_ptr<std::vector<
            std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>
        >> pShadows,
        const int64_t uiRStart,
        const double fAngle
    )
{
    // sort shadows (increasingly) by start coordinate of the match
    std::sort(
            pShadows->begin(),
            pShadows->end(),
            [](std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex> xA,
                std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex> xB)
            {
                /*
                * sort by the interval starts
                * if two intervals start at the same point the larger one shall be treated first
                */
                if(std::get<1>(xA) == std::get<1>(xB))
                    return ( std::get<2>(xA) > std::get<2>(xB) );
                return ( std::get<1>(xA) < std::get<1>(xB) );
            }// lambda
        );// sort function call

    //for(std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>& xTup : *pShadows)
    //{
    //    std::cout << std::get<0>(xTup)->start() << ", " << std::get<0>(xTup)->start_ref() << ", " << std::get<0>(xTup)->size() << std::endl;
    //}

    auto pItervalEnds = std::make_shared<
            std::vector<std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>>
        >();
    pItervalEnds->reserve(pShadows->size());

    nucSeqIndex x = 0;
    //this is the line sweeping part
    for(auto& xTup : *pShadows)
    {
        if(x < std::get<2>(xTup))
        {
            pItervalEnds->push_back(xTup);
            x = std::get<2>(xTup);
        }// if
        else
        {
            assert(! pItervalEnds->empty() );
            double fDistance = deltaDistance(*std::get<0>(xTup), fAngle, uiRStart);;
            //std::cout << "D-";
            nucSeqIndex uiPos = pItervalEnds->size();//uiPos is unsigned!!!
            bool bThisIsCloserToDiagonal = true;
            while(
                    uiPos > 0 && 
                    std::get<2>((*pItervalEnds)[uiPos - 1]) >= std::get<2>(xTup)
                )
            {
                double fDistanceOther = deltaDistance(
                    *std::get<0>((*pItervalEnds)[uiPos - 1]),
                    fAngle,
                    uiRStart
                );
                if(fDistanceOther <= fDistance)
                {
                    bThisIsCloserToDiagonal = false;
                    break;
                }// if
                --uiPos;
            }// while
            if(bThisIsCloserToDiagonal)
            {
                while(
                        ! pItervalEnds->empty() && 
                        std::get<2>(pItervalEnds->back()) >= std::get<2>(xTup)
                    )
                    pItervalEnds->pop_back();
                pItervalEnds->push_back(xTup);
            }// if
            // else do nothing
        }// else
    }// for
    //std::cout << pItervalEnds->size() << std::endl;
    return pItervalEnds;
}// function

std::shared_ptr<Container> LinearLineSweep::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
    const auto& pSoCIn   = std::static_pointer_cast<SoCPriorityQueue>((*vpInput)[0]);
    const auto& pQuery   = std::static_pointer_cast<NucSeq          >((*vpInput)[1]);
#define FILTER_1 ( 0 )
#if FILTER_1
    nucSeqIndex uiAccumulativeSeedLength = 0;
#endif

    unsigned int uiNumTries = 0;
    nucSeqIndex uiLastHarmScore = 0;
    nucSeqIndex uiBestSoCScore = 0;
    unsigned int uiSoCRepeatCounter = 0;


    auto pSoCs = std::make_shared<ContainerVector>(std::make_shared<Seeds>());

    while(!pSoCIn->empty())
    {
        if(++uiNumTries > uiMaxTries)
            break;
        auto pSeedsIn = pSoCIn->pop();

        auto pSeeds = std::make_shared<Seeds>();
        if(pSeedsIn->size() > 1)
        {

            DEBUG(
                pSoCIn->vExtractOrder.push_back(SoCPriorityQueue::blub());
                pSoCIn->vExtractOrder.back().rStartSoC = pSeedsIn->front().start_ref();
                pSoCIn->vExtractOrder.back().rEndSoC = pSeedsIn->front().end_ref();
            )// DEBUG

            assert(!pSeedsIn->empty());
            nucSeqIndex uiCurrSoCScore = 0;
            for(const auto& rSeed : *pSeedsIn)
            {
                uiCurrSoCScore += rSeed.size();
                DEBUG(
                    pSoCIn->vExtractOrder.back().rStartSoC = std::min(
                        pSoCIn->vExtractOrder.back().rStartSoC,
                        rSeed.start_ref()
                    );
                    pSoCIn->vExtractOrder.back().rEndSoC = std::max(
                        pSoCIn->vExtractOrder.back().rEndSoC,
                        rSeed.end_ref()
                    );
                )// DEBUG
            }// for

            // Prof. Kutzners filter:
            // this merely checks weather we actually do have to do the harmonization at all
            if (pQuery->length() > uiSwitchQLen)
            {
                if(uiLastHarmScore > uiCurrSoCScore)
                    continue;
            }// if

            if(uiBestSoCScore * fScoreTolerace > uiCurrSoCScore)
                break;
            uiBestSoCScore = std::max(uiBestSoCScore, uiCurrSoCScore);

            DEBUG(
                pSoCIn->vExtractOrder.back().first = uiCurrSoCScore;
                pSoCIn->vIngroup.push_back(std::make_shared<Seeds>());
            )
#if USE_RANSAC == 1 // switch between ransac line angle + intercept estimation & 45deg median line
            std::vector<double> vX, vY;
            vX.reserve(pSeedsIn->size());
            vY.reserve(pSeedsIn->size());
            for(const auto& rSeed : *pSeedsIn)
            {
                vX.push_back( (double)rSeed.start_ref() + rSeed.size()/2.0);
                vY.push_back( (double)rSeed.start() + rSeed.size()/2.0);

                vX.push_back( (double)rSeed.start_ref());
                vY.push_back( (double)rSeed.start());

                vX.push_back( (double)rSeed.start_ref() + rSeed.size());
                vY.push_back( (double)rSeed.start() + rSeed.size());
            }// for
            /* The Mean Absolute Deviation (MAD) is later required for the threshold t */
            double fMAD = medianAbsoluteDeviation<double>( vY );
            auto xSlopeIntercept = run_ransac(vX, vY, /*pSoCIn->vIngroup.back(),*/ fMAD);

            /*
             * remove outliers
             */
            std::remove_if (
                pSeedsIn->begin(), 
                pSeedsIn->end(), 
                [&]
                (const Seed& rS)
                {
                    return 
                        deltaDistance(
                                rS, 
                                xSlopeIntercept.first, 
                                xSlopeIntercept.second
                            ) > fMAD;
                }
            );
#else
            auto rMedianSeed = (*pSeedsIn)[pSeedsIn->size()/2];
            auto xSlopeIntercept = std::make_pair(
                    0.785398,//forty five degrees
                    (double)rMedianSeed.start_ref() - (double)rMedianSeed.start()
                );
#endif


            DEBUG(
                pSoCIn->vSlopes.push_back(std::tan(xSlopeIntercept.first));
                pSoCIn->vIntercepts.push_back(xSlopeIntercept.second);
            )// DEBUG

            auto pShadows = std::make_shared<
                    std::vector<std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>>
                >();
            pShadows->reserve(pSeedsIn->size());

            // get the left shadows
            for(Seeds::iterator pSeed = pSeedsIn->begin(); pSeed != pSeedsIn->end(); pSeed++)
            {
                pShadows->push_back(std::make_tuple(
                        pSeed,
                        pSeed->start(),
                        pSeed->end_ref()
                    ));
                //std::cout << "(" << pSeed->start() << "," << pSeed->start_ref() << "," << pSeed->size() << ")," << std::endl;
            }//for

            // perform the line sweep algorithm on the left shadows
            auto pShadows2 = linesweep(pShadows, xSlopeIntercept.second, xSlopeIntercept.first);
            pShadows->clear();
            pShadows->reserve(pShadows2->size());

            // get the right shadows
            for(auto &xT : *pShadows2)
                pShadows->push_back(std::make_tuple(
                        std::get<0>(xT),
                        std::get<0>(xT)->start_ref(),
                        std::get<0>(xT)->end()
                    ));

            // perform the line sweep algorithm on the right shadows
            pShadows = linesweep(pShadows, xSlopeIntercept.second, xSlopeIntercept.first);

            pSeeds->reserve(pShadows->size());

            pSeeds->xStats = pSeedsIn->xStats;
            for(auto &xT : *pShadows)
                pSeeds->push_back(*std::get<0>(xT));

            pSeeds->bConsistent = true;

            DEBUG(pSoCIn->vHarmSoCs.push_back(pSeeds);)// DEBUG

            // seeds need to be sorted for the following steps
            std::sort(
                pSeeds->begin(), pSeeds->end(),
                [](const Seed& xA, const Seed& xB)
                {
                    if(xA.start_ref() == xB.start_ref())
                        return xA.start() < xB.start();
                    return xA.start_ref() < xB.start_ref();
                }//lambda
            );//sort function call

            DEBUG(
                pSoCIn->vExtractOrder.back().rStart = pSeeds->front().start_ref();
                pSoCIn->vExtractOrder.back().rEnd = pSeeds->back().end_ref();
            )// DEBUG

            /*
            * sometimes we have one or less seeds remaining after the cupling:
            * in these cases we simply return the center seed (judging by the delta from SoCs)
            * This increases accuracy since the aligner is guided towards the correct
            * position rather than giving it the sound seed or no seed...
            * 
            * Note: returning the longest or leas ambiguous seed does not make sense here;
            * In most cases where this condition triggeres we have seeds that are roughly the same  length
            * and ambiguity... (e.g. 3 seeds of length 17 and two of length 16)
            */
            if(pSeeds->size() <= 1)
            {
                pSeeds->clear();
                pSeeds->push_back( (*pSeedsIn)[pSeedsIn->size()/2]);
            }// if
            assert(!pSeeds->empty());


            /*
             * FILTER:
             * do a gap cost estimation [ O(n) ]
             * this does not improve quality but performance
             * since we might remove some seeds that are too far from another
             *
             * this is much like the backtracking in SW/NW
             * we make sure that we can only have a positive score
             * 
             * Also enforces that no gap is larger than max_gap either on query or reference
             * 
             * Can simply be dis- and enabled, since it's only effect is to delete some seeds.
             */
#if 1
            
            //running score
            int64_t iScore = iMatch * pSeeds->front().size();
            //maximal score
            nucSeqIndex uiMaxScore = iScore;
            //points to the last seed that had a score of 0
            Seeds::iterator pLastStart = pSeeds->begin();
            //outcome
            Seeds::iterator pOptimalStart = pSeeds->begin();
            //outcome
            Seeds::iterator pOptimalEnd = pSeeds->begin();

            /*
            * the goal of this loop is to set pOptimalStart & end correctly
            */
            for(Seeds::iterator pSeed = pSeeds->begin() + 1; pSeed != pSeeds->end(); pSeed++)
            {
                assert(pSeed->start() <= pSeed->end());
                //adjust the score correctly
                iScore += iMatch * pSeed->size();
                /*
                * we need to extract the gap between the seeds in order to do that.
                * we will assume that any gap can be spanned by a single indel
                * while all nucleotides give matches.
                * Therefore we need to get the width x and height y of the rectangular gap
                * number of matches equals min(x, y)
                * gap size equals |x-y|
                */
                nucSeqIndex uiGap = 0;
                if(pSeed->start() > (pSeed-1)->start())
                    uiGap = pSeed->start() - (pSeed-1)->start();
                if(pSeed->start_ref() > (pSeed-1)->start_ref())
                {
                    if(pSeed->start_ref() - (pSeed-1)->start_ref() < uiGap)
                    {
                        uiGap -= pSeed->start_ref() - (pSeed-1)->start_ref();
                        if(optimisticGapEstimation)
                            iScore += iMatch * (pSeed->start_ref() - (pSeed-1)->start_ref());
                    }//if
                    else
                    {
                        if(optimisticGapEstimation)
                            iScore += iMatch * uiGap;
                        uiGap = (pSeed->start_ref() - (pSeed-1)->start_ref()) - uiGap;
                    }//else
                }//if
                uiGap *= iExtend;
                if(uiGap > 0)
                    uiGap += iGap;
                nucSeqIndex uiGapY = pSeed->start() - ( (pSeed-1)->start() + (pSeed-1)->size() );
                if( pSeed->start() < ( (pSeed-1)->start() + (pSeed-1)->size() ) )
                    uiGapY = 0;
                nucSeqIndex uiGapX = pSeed->start_ref() - ((pSeed-1)->start_ref() + (pSeed-1)->size());
                if( pSeed->start_ref() < ( (pSeed-1)->start_ref() + (pSeed-1)->size() ) )
                    uiGapX = 0;
                if( //check for the maximal allowed gap area
                        //uiMaxGapArea == 0 -> disabled
                        ( uiMaxGapArea > 0 && uiGapX > uiMaxGapArea && uiGapY != 0 )
                            ||
                        ( uiMaxGapArea > 0 && uiGapY > uiMaxGapArea && uiGapX != 0  )
                            ||
                        //check for negative score
                        iScore < (int64_t)uiGap
                    )
                {
                    iScore = 0;
                    pLastStart = pSeed;
                }//if
                else
                    iScore -= uiGap;
                if(iScore > (int64_t)uiMaxScore)
                {
                    uiMaxScore = iScore;
                    pOptimalStart = pLastStart;
                    pOptimalEnd = pSeed;
                }//if
            }//for

            /*
            * we need to increment both iterator but check if they extend past the end of pSeeds
            * (check twice because incrementing an iterator pointing to end will
            * result in undefined behaviour)
            *
            * We then simply remove all known suboptimal seeds
            * this is an optimistic estimation so some suboptimal regions might remain
            *
            * NOTE: order is significant!
            * """
            * --- Iterator validity ---
            * Iterators, pointers and references pointing to position (or first) and beyond are
            *             * invalidated, with all iterators, pointers and references to elements before position
            * (or first) are guaranteed to keep referring to the same elements they were referring to
            * before the call.
            * """
            */
            if(pOptimalEnd != pSeeds->end())
                if(++pOptimalEnd != pSeeds->end())
                    pSeeds->erase(pOptimalEnd, pSeeds->end());
            if(pOptimalStart != pSeeds->end())
                pSeeds->erase(pSeeds->begin(), pOptimalStart);
#endif

        }// if
        else // pSeedsIn contains merely one seed
        {
            pSeeds->push_back(pSeedsIn->front());
            DEBUG(
                pSoCIn->vExtractOrder.push_back(SoCPriorityQueue::blub());
                pSoCIn->vExtractOrder.back().first = pSeedsIn->front().size();
                pSoCIn->vSlopes.push_back(0);
                pSoCIn->vIntercepts.push_back(0);
                pSoCIn->vHarmSoCs.push_back(std::make_shared<Seeds>(pSeeds));
            )// DEBUG
        }

        /*
         * end of FILTER
         */
        nucSeqIndex uiCurrHarmScore = 0;
        for(const auto& rSeed : *pSeeds)
            uiCurrHarmScore += rSeed.size();

        if(uiCurrHarmScore < 18 )
            continue;
        if(uiCurrHarmScore < pQuery->length() * 0.002 )
            continue;

        DEBUG(
            std::vector<bool> vQCoverage(pQuery->length(), false);
            pSoCIn->vExtractOrder.back().qCoverage = 0;
            for(const auto& rSeed : *pSeeds)
            {
                for(auto uiX = rSeed.start(); uiX < rSeed.end(); uiX++)
                    vQCoverage[uiX] = true;
            }//for
            pSoCIn->vExtractOrder.back().second = uiCurrHarmScore;

            for(bool b : vQCoverage)
                if(b)
                    pSoCIn->vExtractOrder.back().qCoverage++;
        )// DEBUG
#if 1
        // Prof. Kutzners filter: (is this equivalent to just always looking at the best SoC..?)
        if (pQuery->length() > uiSwitchQLen)
        {
            if(uiLastHarmScore > uiCurrHarmScore)
                continue;
        }// if
        else
        {
            if(
                !(uiCurrHarmScore + (pQuery->length()*fScoreDiffTolerance) >= uiLastHarmScore &&
                uiCurrHarmScore - (pQuery->length()*fScoreDiffTolerance) <= uiLastHarmScore)
                )
                uiSoCRepeatCounter = 0;
            else if(++uiSoCRepeatCounter >= uiMaxEqualScoreLookahead)
            {
                uiSoCRepeatCounter -= 1; // cause we haven't actually pushed the current soc yet...
                break;
            }// else
        }// else
#endif
        uiLastHarmScore = uiCurrHarmScore;

        pSoCs->push_back(pSeeds);

        //FILTER
#if FILTER_1

        nucSeqIndex uiAccLen = pSeeds->getScore();
        if (uiAccumulativeSeedLength > uiAccLen )
            continue;
        uiAccumulativeSeedLength = std::max(uiAccLen, uiAccumulativeSeedLength);
#endif
        //FILTER END

    }//while

    for(unsigned int ui = 0; ui < uiSoCRepeatCounter && pSoCs->size() > 1; ui++)
        pSoCs->pop_back();
    
    return pSoCs;
}//function

#ifdef WITH_PYTHON
void exportLinesweep()
{
    //export the LineSweepContainer class
    boost::python::class_<
            LinearLineSweep, 
            boost::python::bases<Module>,
            std::shared_ptr<LinearLineSweep>
        >("LinearLineSweep")
        .def_readwrite("optimistic_gap_estimation", &LinearLineSweep::optimisticGapEstimation)
        .def_readwrite("min_coverage", &LinearLineSweep::fMinimalQueryCoverage)
        .def_readwrite("tolerance", &LinearLineSweep::fScoreTolerace)
        .def_readwrite("max_tries", &LinearLineSweep::uiMaxTries)
        .def_readwrite("equal_score_lookahead", &LinearLineSweep::uiMaxEqualScoreLookahead)
        .def_readwrite("diff_tolerance", &LinearLineSweep::fScoreDiffTolerance)
        .def_readwrite("switch_q_len", &LinearLineSweep::uiSwitchQLen)
    ;
    boost::python::implicitly_convertible< 
        std::shared_ptr<LinearLineSweep>,
        std::shared_ptr<Module> 
    >();
}//function
#endif