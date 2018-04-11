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
extern float fRelativePadding;


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
            while(!pItervalEnds->empty() && std::get<2>(pItervalEnds->back()) >= std::get<2>(xTup))
                pItervalEnds->pop_back();
    }//for
    return pItervalEnds;
}//function

nucSeqIndex TempBackend::needlemanWunsch(
        std::shared_ptr<NucSeq> pQuery, 
        std::shared_ptr<NucSeq> pRef,
        nucSeqIndex fromQuery,
        nucSeqIndex toQuery,
        nucSeqIndex fromRef,
        nucSeqIndex toRef,
        std::shared_ptr<Alignment> pAlignment,
        bool bNoGapAtBeginning,
        bool bNoGapAtEnd
        DEBUG_PARAM(bool bPrintMatrix)
    )
{
    // in these cases parasail does not offer the wanted alignment approach
    // so we use our own slow implementation
    if(bNoGapAtBeginning || bNoGapAtEnd)
        return naiveNeedlemanWunsch(
            pQuery, pRef,
            fromQuery, toQuery,
            fromRef, toRef,
            pAlignment,
            bNoGapAtBeginning, bNoGapAtEnd
            DEBUG_PARAM(bPrintMatrix)
        );
    // in all other cases we use parasail
    // do some checking for empty sequences though since parasail does not offer that
    if(toRef <= fromRef)
        if(toQuery <= fromQuery)
            return 0;
    if(toQuery <= fromQuery)
    {
        pAlignment->append(MatchType::deletion, toRef-fromRef);
        return 0;
    }//if
    if(toRef <= fromRef)
    {
        pAlignment->append(MatchType::insertion, toQuery-fromQuery);
        return 0;
    }//if

    // okay if we reached here we actually have to align something

    /*
     * do the NW alignment
     */

    // Note: parasail does not follow the usual theme where for opening a gap 
   ParsailResultWrapper pResult(parasail_nw_trace_scan_16(
        (const char*)pQuery->pGetSequenceRef() + fromQuery, toQuery - fromQuery,
        (const char*)pRef->pGetSequenceRef() + fromRef, toRef - fromRef,
        iGap + iExtend, iExtend, &matrix
    ));

    //get the cigar
    ParsailCigarWrapper pCigar(parasail_result_get_cigar(
        pResult.get(),
        (const char*)pQuery->pGetSequenceRef() + fromQuery, toQuery - fromQuery,
        (const char*)pRef->pGetSequenceRef() + fromRef, toRef - fromRef,
        &matrix
    ));

    DEBUG_2(
        std::cout << "cigar length: " << pCigar->len << std::endl;
        std::cout << pCigar->beg_query << ", " << pCigar->beg_ref << std::endl;
    )

    //Due to using NW this should hold:
    assert(pCigar->beg_query == 0);
    assert(pCigar->beg_ref == 0);

    /*
     * Decode the cigar
     */

    nucSeqIndex uiQPos = 0;
    nucSeqIndex uiRPos = 0;
    //decode the cigar
    for(int i = 0; i < pCigar->len; i++)
    {
        char c = parasail_cigar_decode_op(pCigar->seq[i]);
        uint32_t uiLen = parasail_cigar_decode_len(pCigar->seq[i]);
        DEBUG_2(std::cout << c << " x" << uiLen << std::endl;)
        switch (c)
        {
            case '=':
                pAlignment->append(MatchType::match, uiLen);
                uiQPos += uiLen;
                uiRPos += uiLen;
                break;
            case 'I':
                pAlignment->append(MatchType::insertion, uiLen);
                uiQPos += uiLen;
                break;
            case 'D':
                pAlignment->append(MatchType::deletion, uiLen);
                uiRPos += uiLen;
                break;
            case 'X':
                pAlignment->append(MatchType::missmatch, uiLen);
                uiQPos += uiLen;
                uiRPos += uiLen;
                break;
            default:
                // there are different CIGAR symbols allowed in the SAM format
                // but parasail should never generate any of them
                assert(false);
                break;
        }//switch
    }//for
    
    //Due to using NW this should hold:
    assert(uiQPos == toQuery - fromQuery);
    assert(uiRPos == toRef - fromRef);
    
    DEBUG_2(
        std::cout << pQuery->length() - uiQPos << ", " << pRef->length() - uiRPos << std::endl;
    )

    //since we use the naive implementation for semi-global alignments we can always return 0 here
    return 0;
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

    while(!pSoCIn->empty())
    {
        //@note from here on it is the original linesweep module
        auto pSeedsIn = pSoCIn->pop();

        auto pShadows = std::make_shared<
                std::vector<std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>>
            >();

        // get the left shadows
        for(Seeds::iterator pSeed = pSeedsIn->begin(); pSeed != pSeedsIn->end(); pSeed++)
            pShadows->push_back(std::make_tuple(
                    pSeed,
                    pSeed->start(),
                    pSeed->end_ref()
                ));

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
        * FILTER:
        * do a gap cost estimation [O(n)]
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
            assert(pSeed->size() <= 10000);
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

        /*
        * One more thing: 
        * sometimes we have one or less seeds remaining after the cupling:
        * in these cases we simply return the center seed (judging by the delta from SoCs)
        * This increases accuracy since the aligner is guided towards the correct
        * position rather than giving it the sound seed or no seed...
        * 
        * Note: returning the longest or leas ambiguous seed does not make sense here;
        * In most cases where this condition triggeres we have seeds that are roughly the same length
        * and ambiguity... (e.g. 3 seeds of length 17 and two of length 16)
        */
        if(pSeeds->size() <= 1)
            if(!pSeedsIn->empty())
                pSeeds->push_back( (*pSeedsIn)[pSeedsIn->size()/2]);

        assert(!pSeeds->empty());

        //@note from here on it is the original NW module

        //making a lambda function so that i do not have to edit the original code
        auto fNWModule = [&]()
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

            /*
            * Here we decide weather to actually perform NW in the gaps between seeds or 
            * if we use SW to align the entire thing.
            */
            if(endQuery - beginQuery < pQuery->length() * fMinimalQueryCoverage)
            {
                DEBUG_2(std::cout << "computing SW for entire area" << std::endl;)
                //figure out the correct reference location
                nucSeqIndex refStart = beginRef;
                if(refStart >= pQuery->length() * fRelativePadding)
                    refStart -=  pQuery->length() * fRelativePadding;
                else
                    refStart = 0;
                nucSeqIndex refEnd = endRef;
                assert(endQuery <= pQuery->length());
                refEnd += pQuery->length()  * fRelativePadding;

                assert(refEnd > refStart);
                nucSeqIndex refWidth = refEnd - refStart;
                if(refStart + refWidth >= pRefPack->uiUnpackedSizeForwardPlusReverse())
                    refWidth = pRefPack->uiUnpackedSizeForwardPlusReverse() - refStart - 1;
                assert(refWidth > 0);
                if(pRefPack->bridgingSubsection(refStart, refWidth))
                {
                    DEBUG_2(std::cout << "Un-bridging from " << refStart << ", " << refWidth;)
                    pRefPack->unBridgeSubsection(refStart, refWidth);
                    DEBUG_2(std::cout << " to: " << refStart << ", " << refWidth << std::endl;)
                }
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

            // here we have enough query coverage to attemt to fill in the gaps merely
            DEBUG_2(std::cout << "filling in gaps" << std::endl;)

            std::shared_ptr<Alignment> pRet;

            if(!bLocal)
            {
                beginRef -= (nucSeqIndex)(  pQuery->length() * fRelativePadding );
                if(beginRef > endRef)//check for underflow
                    beginRef = 0;
                endRef += (nucSeqIndex)( pQuery->length() * fRelativePadding );
                if(beginRef > endRef)//check for overflow
                    endRef = pRefPack->uiUnpackedSizeForwardPlusReverse()-1;
                endQuery = pQuery->length();
                beginQuery = 0;
            }//if
            assert(endQuery <= pQuery->length());
            pRet = std::shared_ptr<Alignment>(
                new Alignment(beginRef, endRef, beginQuery, endQuery)
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
            }

            //create the actual alignment

            if(!bLocal)
            {
                pRet->uiBeginOnRef += needlemanWunsch(
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
                //skip the first seed
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
                    needlemanWunsch(
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
                pRet->uiEndOnRef -= needlemanWunsch(
                    pQuery,
                    pRef,
                    endOfLastSeedQuery,
                    endQuery-1,
                    endOfLastSeedReference,
                    endRef-beginRef-1,
                    pRet,
                    false,
                    true
                );
                pRet->removeDangeling();
            }//else
            return pRet;

        };// lambda function

        //get the NW alignment
        std::shared_ptr<Alignment> pAlignment = fNWModule();


        //@note from here on is the new code

        // check weather the newly discovered alignment is the best one
        if( pBestAlignment->score() < pAlignment->score() )
            pBestAlignment = pAlignment;
        pAlignmentsOut->push_back(pAlignment);

        if(pAlignmentsOut->size() >= uiMaxTries)
            break;

        if( pBestAlignment->score() * fScoreTolerace > pAlignment->score() )
            //we found an alignment where we a reasonably certain that it is the best one
            break;
    }//while

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
    ;
    boost::python::implicitly_convertible< 
        std::shared_ptr<TempBackend>,
        std::shared_ptr<Module> 
    >();
}//function

