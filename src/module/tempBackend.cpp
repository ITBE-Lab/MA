/** 
 * @file tempBackend.cpp
 * @author Markus Schmidt
 */
#include "module/tempBackend.h"
#include "module/needlemanWunsch.h"
#include "container/soc.h"
#include "container/pack.h"
using namespace libMA;

extern int iGap;
extern int iExtend;
extern int iMatch;
extern int iMissMatch;
extern nucSeqIndex uiMaxGapArea;
/*
 * @todo: fix this problem
 * 
 * arghh this is really ugly...
 * At the moment there is only one sequence that sets all NW and SW parameters correctly:
 * 1) set the parameters.
 * 2) create a new NW module.
 * 3) Then create all other modules that you want to use.....
 */
//the match missmatch matrix
extern parasail_matrix_t matrix;
extern std::vector<int> vMatrixContent;
extern const int parasail_custom_map[];
extern nucSeqIndex uiPadding;


ContainerVector TempBackend::getInputType() const
{
    return ContainerVector{
            std::make_shared<SoCPriorityQueue>(),
            std::make_shared<NucSeq>(),
            std::make_shared<Pack>()
        };
}//function

std::shared_ptr<Container> TempBackend::getOutputType() const
{
    return std::make_shared<ContainerVector>( std::make_shared<Alignment>() );
}//function

std::shared_ptr<std::vector<std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>>>
    TempBackend::linesweep(
        std::shared_ptr<std::vector<
            std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>
        >> pShadows
    )
{
    //sort shadows (increasingly) by start coordinate of the match
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
                    return std::get<2>(xA) > std::get<2>(xB);
                return std::get<1>(xA) < std::get<1>(xB);
            }//lambda
        );//sort function call

    //for(std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>& xTup : *pShadows)
    //{
    //    std::cout << std::get<0>(xTup)->start() << ", " << std::get<0>(xTup)->start_ref() << ", " << std::get<0>(xTup)->size() << std::endl;
    //}

    auto pItervalEnds = std::make_shared<
            std::vector<std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>>
        >();

    nucSeqIndex x = 0;
    //this is the line sweeping part
    for(std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>& xTup : *pShadows)
    {
        if(x < std::get<2>(xTup))
        {
            pItervalEnds->push_back(xTup);
            x = std::get<2>(xTup);
        }//if
        /*
         * sometimes a seeding algorithm may deliver equal seeds
         * in that case these seeds are obviously not contradicting
         *
         * This is not mentioned in the paper since we talk of seed SETS there.
         * A set cannot have the same element twice,
         * however the vector we use in practice can...
         *
         * Anyways we just need to add one simple check for equality here
         * operator!= is not implemented so we rely on ! of ==.
         */
        else if(
                ! pItervalEnds->empty() 
                && ! (*std::get<0>(xTup) == *std::get<0>(pItervalEnds->back())) 
            )
        {
            //std::cout << "D-";
            while(!pItervalEnds->empty() && std::get<2>(pItervalEnds->back()) >= std::get<2>(xTup))
            {
                //std::cout << "d-";
                pItervalEnds->pop_back();
            }/// while
        }//else if
    }//for
    //std::cout << std::endl;
    return pItervalEnds;
}//function

std::shared_ptr<Container> TempBackend::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
    const auto& pSoCIn   = std::static_pointer_cast<SoCPriorityQueue>((*vpInput)[0]);
    const auto& pQuery   = std::static_pointer_cast<NucSeq          >((*vpInput)[1]);
    const auto& pRefPack = std::static_pointer_cast<Pack            >((*vpInput)[2]);

    auto pAlignmentsOut = std::make_shared<ContainerVector>(std::make_shared<Alignment>());
    auto pBestAlignment = std::make_shared<Alignment>();

#define FILTER_1 ( 0 )
#if FILTER_1
    nucSeqIndex uiAccumulativeSeedLength = 0;
#endif

    unsigned int uiNumTries = 0;
    nucSeqIndex uiLastHarmScore = 0;
    nucSeqIndex uiBestSoCScore = 0;
    unsigned int uiSoCRepeatCounter = 0;

    std::vector<std::shared_ptr<Seeds>> vSoCs;

    while(!pSoCIn->empty())
    {
        if(++uiNumTries > uiMaxTries)
            break;
        //@note from here on it is the original linesweep module
        auto pSeedsIn = pSoCIn->pop();
        assert(!pSeedsIn->empty());
        nucSeqIndex uiCurrSoCScore = 0;
        for(const auto& rSeed : *pSeedsIn)
            uiCurrSoCScore += rSeed.size();

        DEBUG(
            pSoCIn->vExtractOrder.push_back(SoCPriorityQueue::blub());
            pSoCIn->vExtractOrder.back().first = uiCurrSoCScore;
        )// DEBUG

        auto pShadows = std::make_shared<
                std::vector<std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>>
            >();

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
        auto pShadows2 = linesweep(pShadows);
        pShadows->clear();

        // get the right shadows
        for(auto &xT : *pShadows2)
            pShadows->push_back(std::make_tuple(
                    std::get<0>(xT),
                    std::get<0>(xT)->start_ref(),
                    std::get<0>(xT)->end()
                ));

        // perform the line sweep algorithm on the right shadows
        pShadows = linesweep(pShadows);

        auto pSeeds = std::make_shared<Seeds>();
        pSeeds->xStats = pSeedsIn->xStats;
        for(auto &xT : *pShadows)
            pSeeds->push_back(*std::get<0>(xT));

        pSeeds->bConsistent = true;

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
         */
        //query position of last encountered seed
        nucSeqIndex uiLastQ = pSeeds->front().start();
        //reference position of last encountered seed
        nucSeqIndex uiLastR = pSeeds->front().start_ref();
        //running score
        unsigned long uiScore = 0;
        //maximal score
        nucSeqIndex uiMaxScore = 0;
        //points to the last seed that had a score of 0
        Seeds::iterator pLastStart = pSeeds->begin();
        //outcome
        Seeds::iterator pOptimalStart = pSeeds->begin();
        //outcome
        Seeds::iterator pOptimalEnd = pSeeds->end();

        /*
         * the goal of this loop is to set pOptimalStart & end correctly
         */
        for(Seeds::iterator pSeed = pSeeds->begin(); pSeed != pSeeds->end(); pSeed++)
        {
            assert(pSeed->start() <= pSeed->end());
            //adjust the score correctly
            uiScore += iMatch * pSeed->size();
            /*
             * we need to extract the gap between the seeds in order to do that.
             * we will assume that any gap can be spanned by a single indel
             * while all nucleotides give matches.
             * Therefore we need to get the width x and height y of the rectangular gap
             * number of matches equals min(x, y)
             * gap size equals |x-y|
             */
            nucSeqIndex uiGap = 0;
            if(pSeed->start() > uiLastQ)
                uiGap = pSeed->start() - uiLastQ;
            if(pSeed->start_ref() > uiLastR)
            {
                if(pSeed->start_ref() - uiLastR < uiGap)
                {
                    uiGap -= pSeed->start_ref() - uiLastR;
                    if(optimisticGapEstimation)
                        uiScore += iMatch * (pSeed->start_ref() - uiLastR);
                }//if
                else
                {
                    if(optimisticGapEstimation)
                        uiScore += iMatch * uiGap;
                    uiGap = (pSeed->start_ref() - uiLastR) - uiGap;
                }//else
            }//if
            uiGap *= iExtend;
            if(uiGap > 0)
                uiGap += iGap;
            if( //check for the maximal allowed gap area
                (uiMaxGapArea > 0 && //0 == disabled
                    (pSeed->start() >= uiLastQ ? pSeed->start() - uiLastQ : 1) *
                    (pSeed->start_ref() >= uiLastR ? pSeed->start_ref() - uiLastR : 1)
                        > uiMaxGapArea) 
                    || 
                //check for negative score
                uiScore < uiGap)
            {
                uiScore = 0;
                pLastStart = pSeed;
            }//if
            else
                uiScore -= uiGap;
            if(uiScore > uiMaxScore)
            {
                uiMaxScore = uiScore;
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
         * invalidated, with all iterators, pointers and references to elements before position
         * (or first) are guaranteed to keep referring to the same elements they were referring to
         * before the call.
         * """
         */
        if(pOptimalEnd != pSeeds->end())
            if(++pOptimalEnd != pSeeds->end())
                pSeeds->erase(pOptimalEnd, pSeeds->end());
        if(pOptimalStart != pSeeds->end())
            if(++pOptimalStart != pSeeds->end())
                pSeeds->erase(pSeeds->begin(), pOptimalStart);
        /*
         * end of FILTER
         */
        nucSeqIndex uiCurrHarmScore = 0;
        for(const auto& rSeed : *pSeeds)
            uiCurrHarmScore += rSeed.size();

        // Prof. Kutzners killer filter:
        if (pQuery->length() > 500 && uiLastHarmScore > uiCurrHarmScore)
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

        if(
            !(uiCurrHarmScore + (pQuery->length()*fScoreDiffTolerance) >= uiLastHarmScore &&
              uiCurrHarmScore - (pQuery->length()*fScoreDiffTolerance) <= uiLastHarmScore)
            )
            uiSoCRepeatCounter = 0;
        else if(++uiSoCRepeatCounter >= uiMaxEqualScoreLookahead)
        {
            uiSoCRepeatCounter -= 1; // cause we havent actually pushed the current soc yet...
            break;
        }//else

        uiLastHarmScore = uiCurrHarmScore;

        if(uiBestSoCScore * fScoreTolerace > uiCurrSoCScore)
            break;

        uiBestSoCScore = std::max(uiBestSoCScore, uiCurrSoCScore);

        vSoCs.push_back(pSeeds);

        //FILTER
#if FILTER_1

        nucSeqIndex uiAccLen = pSeeds->getScore();
        if (uiAccumulativeSeedLength > uiAccLen )
            continue;
        uiAccumulativeSeedLength = std::max(uiAccLen, uiAccumulativeSeedLength);
#endif
        //FILTER END

    }//while

    for(unsigned int ui = 0; ui < uiSoCRepeatCounter && !vSoCs.empty(); ui++)
        vSoCs.pop_back();

    DEBUG(
        unsigned int uiDCounter = 0;
    )//DEBUG

    for (auto pSeeds : vSoCs)
    {
        //@note from here on it is the original NW module

        //making a lambda function so that i do not have to edit the original code
        const auto fNWModule = [&]()
        {
            if(pSeeds == nullptr)
                return std::shared_ptr<Alignment>(new Alignment());

            //no seeds => no spot found at all...
            if(pSeeds->empty())
            {
                std::shared_ptr<Alignment> pRet(new Alignment());
                pRet->xStats = pSeeds->xStats;
                pRet->xStats.sName = pQuery->sName;
                return pRet;
            }//if

            DEBUG_2(
                std::cout << "seedlist: (start_ref, end_ref; start_query, end_query)" << std::endl;
                for(Seed& rSeed : *pSeeds)
                {
                    std::cout << rSeed.start_ref() << ", " << rSeed.end_ref() << "; "
                        << rSeed.start() << ", " << rSeed.end() << std::endl;
                }//for
            )// DEBUG

            // Determine the query and reverence coverage of the seeds

            nucSeqIndex beginRef = pSeeds->front().start_ref();
            nucSeqIndex endRef = pSeeds->back().end_ref();
            //seeds are sorted by ther startpos so we 
            //actually need to check all seeds to get the proper end
            nucSeqIndex endQuery = pSeeds->back().end();
            nucSeqIndex beginQuery = pSeeds->front().start();
            for (auto xSeed : *pSeeds)
            {
                if(endRef < xSeed.end_ref())
                    endRef = xSeed.end_ref();
                if(beginRef > xSeed.start_ref())
                    beginRef = xSeed.start_ref();
                if(endQuery < xSeed.end())
                    endQuery = xSeed.end();
                if(beginQuery > xSeed.end())
                    beginQuery = xSeed.start();
                assert(xSeed.start() <= xSeed.end());
            }//for
            DEBUG_2(
                std::cout << beginRef << ", " << endRef << "; " << beginQuery << ", " << endQuery << std::endl;
            )// DEEBUG
#if 0
            /*
            * Here we decide weather to actually perform NW in the gaps between seeds or 
            * if we use SW to align the entire thing.
            */
            if(endQuery - beginQuery < pQuery->length() * fMinimalQueryCoverage)
            {
                DEBUG_2(std::cout << "computing SW for entire area" << std::endl;)
                //figure out the correct reference location
                nucSeqIndex refStart = beginRef;
                if(refStart >= uiPadding + beginQuery)
                    refStart -= uiPadding + beginQuery;
                else
                    refStart = 0;
                nucSeqIndex refEnd = endRef;
                assert(endQuery <= pQuery->length());
                refEnd += uiPadding + (pQuery->length() - endQuery);

                assert(refEnd >= refStart);
                nucSeqIndex refWidth = refEnd - refStart + 1;
                if(refStart + refWidth >= pRefPack->uiUnpackedSizeForwardPlusReverse())
                    refWidth = pRefPack->uiUnpackedSizeForwardPlusReverse() - refStart - 1;
                assert(refWidth > 0);
                if(pRefPack->bridgingSubsection(refStart, refWidth))
                {
                    DEBUG_2(std::cout << "Un-bridging from " << refStart << ", " << refWidth;)
                    pRefPack->unBridgeSubsection(refStart, refWidth);
                    DEBUG_2(std::cout << " to: " << refStart << ", " << refWidth << std::endl;)
                }

                //std::cout << refStart << " to " << refStart + refWidth << std::endl;

                //extract the reference
                std::shared_ptr<NucSeq> pRef = pRefPack->vExtract(
                    refStart,
                    refStart + refWidth
                );
                //compute the SW alignment
                auto pRet = smithWaterman(pQuery, pRef, refStart);
                //copy the stats
                pRet->xStats = pSeeds->xStats;
                pRet->xStats.sName = pQuery->sName;
                // return the SW alignment
                return pRet;
            }//if
#endif

            // here we have enough query coverage to attemt to fill in the gaps merely
            DEBUG_2(std::cout << "filling in gaps" << std::endl;)

            std::shared_ptr<Alignment> pRet;

            if(!bLocal)
            {
                beginRef -= uiPadding + beginQuery;
                if(beginRef > endRef)//check for underflow
                    beginRef = 0;
                endRef += uiPadding + (pQuery->length() - endQuery);
                if(beginRef > endRef)//check for overflow
                    endRef = pRefPack->uiUnpackedSizeForwardPlusReverse();
                endQuery = pQuery->length();
                beginQuery = 0;
            }//if
            assert(endQuery <= pQuery->length());
            pRet = std::shared_ptr<Alignment>(
                new Alignment(beginRef, beginQuery)
            );

            //save the strip of consideration stats in the alignment
            pRet->xStats = pSeeds->xStats;
            pRet->xStats.sName = pQuery->sName;

            DEBUG_2(
                std::cout << beginRef << " " << endRef << std::endl;
            )
            std::shared_ptr<NucSeq> pRef;
            try
            {
                pRef = pRefPack->vExtract(beginRef, endRef);
            } catch(std::runtime_error e)
            {
                std::shared_ptr<Alignment> pRet(new Alignment());
                pRet->xStats = pSeeds->xStats;
                pRet->xStats.sName = pQuery->sName;
                return pRet;
            }// catch

            //create the actual alignment

            if(!bLocal)
            {
                dynPrg(
                    pQuery,
                    pRef,
                    0,
                    pSeeds->front().start(),
                    0,
                    pSeeds->front().start_ref() - beginRef,
                    pRet,
                    true,
                    false
                );
            }//else

            nucSeqIndex endOfLastSeedQuery = pSeeds->front().end();
            nucSeqIndex endOfLastSeedReference = pSeeds->front().end_ref() - beginRef;

            pRet->append(MatchType::seed, pSeeds->front().size());
            bool bSkip = true;
            for(Seed& rSeed : *pSeeds)
            {
                // skip the first seed
                // we do this since the seed has already been appended before the loop
                // this makes the loop structure easier since this way 
                // we can always first compute the NW and then append a seed
                if(bSkip)
                {
                    bSkip = false;
                    continue;
                }//if
                nucSeqIndex ovQ = endOfLastSeedQuery - rSeed.start();
                if(rSeed.start() > endOfLastSeedQuery)
                    ovQ = 0;
                nucSeqIndex ovR = endOfLastSeedReference - (rSeed.start_ref() - beginRef);
                if(rSeed.start_ref() > endOfLastSeedReference + beginRef)
                    ovR = 0;
                nucSeqIndex len = rSeed.size();
                nucSeqIndex overlap = std::max(ovQ, ovR);
                DEBUG_2(
                    std::cout << "overlap: " << overlap << std::endl;
                )//DEBUG
                if(len > overlap)
                {
                    dynPrg(
                            pQuery,
                            pRef,
                            endOfLastSeedQuery,
                            rSeed.start(),
                            endOfLastSeedReference,
                            rSeed.start_ref() - beginRef,
                            pRet,
                            false,
                            false
                        );
                    if(ovQ > ovR)
                        pRet->append(MatchType::deletion, ovQ - ovR);
                    DEBUG_2(
                        for(nucSeqIndex i = ovR; i < ovQ; i++)
                            std::cout << "d";
                    )
                    if(ovR > ovQ)
                        pRet->append(MatchType::insertion, ovR - ovQ);
                    DEBUG_2(
                        for(nucSeqIndex i = ovQ; i < ovR; i++)
                            std::cout << "i";
                        std::cout << std::endl;
                    )//DEBUG
                    pRet->append(MatchType::seed, len - overlap);
                    DEBUG_2(
                        std::cout << len - overlap << std::endl;
                    )//DEBUG_2
                    DEBUG_2(
                        for(nucSeqIndex i = overlap; i < len; i++)
                            std::cout << pQuery->charAt(i + rSeed.start());
                        std::cout << std::endl;
                        for(nucSeqIndex i = overlap; i < len; i++)
                            std::cout << pRef->charAt(i + rSeed.start_ref() - beginRef);
                        std::cout << std::endl;
                    )//DEBUG
                    DEBUG_2(
                        for(nucSeqIndex i = 0; i < len - overlap; i++)
                            std::cout << "m";
                    )//DEBUG_2
                    if(rSeed.end() > endOfLastSeedQuery)
                        endOfLastSeedQuery = rSeed.end();
                    if(rSeed.end_ref() > endOfLastSeedReference + beginRef)
                        endOfLastSeedReference = rSeed.end_ref() - beginRef;
                }//if
            }//for

            if(bLocal)
                assert(std::get<0>(pRet->data.front()) == MatchType::seed);
            assert(std::get<0>(pRet->data.back()) == MatchType::seed);

            DEBUG_2(
                std::cout << std::endl;
            )
            if(bLocal)
                pRet->makeLocal();
            else
            {
                dynPrg(
                    pQuery,
                    pRef,
                    endOfLastSeedQuery,
                    endQuery,
                    endOfLastSeedReference,
                    endRef-beginRef,
                    pRet,
                    false,
                    true
                );
                pRet->removeDangeling();
            }//else
            return pRet;

        };// lambda function

        //this is the NW module
        std::shared_ptr<Alignment> pAlignment = fNWModule();


        //@note from here on is the new code

        // check weather the newly discovered alignment is the best one
        if( pBestAlignment->score() < pAlignment->score() )
            pBestAlignment = pAlignment;
        pAlignmentsOut->push_back(pAlignment);

        DEBUG(
            pSoCIn->vExtractOrder[uiDCounter++].third = pAlignment->score();
        )// DEBUG

        if(pAlignmentsOut->size() >= uiMaxTries)
            break;

    }// for

    //sort alignments ascending
    std::sort(
        pAlignmentsOut->begin(), pAlignmentsOut->end(),
        []
        (std::shared_ptr<Container> a, std::shared_ptr<Container> b)
        {
            return a->larger(b);
        }//lambda
    );//sort function call
    assert(pAlignmentsOut->size() <= 1 || !pAlignmentsOut->back()->larger(pAlignmentsOut->front()));

    return pAlignmentsOut;
}//function

void exportTempBackend()
{
    //export the LineSweepContainer class
    boost::python::class_<
            TempBackend, 
            boost::python::bases<Module>,
            std::shared_ptr<TempBackend>
        >("TempBackend")
            .def_readwrite("min_coverage", &TempBackend::fMinimalQueryCoverage)
            .def_readwrite("tolerance", &TempBackend::fScoreTolerace)
            .def_readwrite("local", &TempBackend::bLocal)
            .def_readwrite("max_tries", &TempBackend::uiMaxTries)
            .def_readwrite("equal_score_lookahead", &TempBackend::uiMaxEqualScoreLookahead)
            .def_readwrite("diff_tolerance", &TempBackend::fScoreDiffTolerance)
    ;
    boost::python::implicitly_convertible< 
        std::shared_ptr<TempBackend>,
        std::shared_ptr<Module> 
    >();
}//function

