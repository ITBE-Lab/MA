/**
 * @file harmonization.cpp
 * @author Markus Schmidt
 */
#include "ma/module/harmonization.h"
#include "ma/module/stripOfConsideration.h"
#include "ms/util/pybind11.h"
#if USE_RANSAC == 1
#include "ma/sample_consensus/test_ransac.h"
#endif
using namespace libMA;
using namespace libMS;

std::shared_ptr<Seeds> Harmonization::applyFilters( std::shared_ptr<Seeds>& pIn ) const
{
    auto pRet = std::make_shared<Seeds>( );
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
     */
    if( bDoGapCostEstimationCutting )
    {
        // running score
        int64_t iScore = pGlobalParams->iMatch->get( ) * pIn->front( ).size( );
        // maximal score
        nucSeqIndex uiMaxScore = iScore;
        // points to the last seed that had a score of 0
        Seeds::iterator pLastStart = pIn->begin( );
        // outcome
        Seeds::iterator pOptimalStart = pIn->begin( );
        // outcome
        Seeds::iterator pOptimalEnd = pIn->begin( );

        /*
         * the goal of this loop is to set pOptimalStart & end correctly
         */
        for( Seeds::iterator pSeed = pIn->begin( ) + 1; pSeed != pIn->end( ); pSeed++ )
        {
            assert( pSeed->start( ) <= pSeed->end( ) );
            // adjust the score correctly
            iScore += pGlobalParams->iMatch->get( ) * pSeed->size( );
            /*
             * we need to extract the gap between the seeds in order to do that.
             * we will assume that any gap can be spanned by a single indel
             * while all nucleotides give matches.
             * Therefore we need to get the width x and height y of the rectangular gap
             * number of matches equals min(x, y)
             * gap size equals |x-y|
             */
            nucSeqIndex uiGap = 0;
            if( pSeed->start( ) > ( pSeed - 1 )->start( ) )
                uiGap = pSeed->start( ) - ( pSeed - 1 )->start( );
            if( pSeed->start_ref( ) > ( pSeed - 1 )->start_ref( ) )
            {
                if( pSeed->start_ref( ) - ( pSeed - 1 )->start_ref( ) < uiGap )
                {
                    uiGap -= pSeed->start_ref( ) - ( pSeed - 1 )->start_ref( );
                    if( optimisticGapEstimation )
                        iScore += pGlobalParams->iMatch->get( ) * ( pSeed->start_ref( ) - ( pSeed - 1 )->start_ref( ) );
                } // if
                else
                {
                    if( optimisticGapEstimation )
                        iScore += pGlobalParams->iMatch->get( ) * uiGap;
                    uiGap = ( pSeed->start_ref( ) - ( pSeed - 1 )->start_ref( ) ) - uiGap;
                } // else
            } // if
            uiGap *= pGlobalParams->iExtend->get( );
            if( uiGap > 0 )
                uiGap += pGlobalParams->iGap->get( );
            if( uiGap > (nucSeqIndex)pGlobalParams->uiSVPenalty->get( ) && pGlobalParams->uiSVPenalty->get( ) != 0 )
                uiGap = (nucSeqIndex)pGlobalParams->uiSVPenalty->get( );
            if( // check for the maximal allowed gap area
                // check for negative score
                iScore < (int64_t)uiGap )
            {
                iScore = 0;
                pLastStart = pSeed;
            } // if
            else
                iScore -= uiGap;
            if( iScore > (int64_t)uiMaxScore )
            {
                uiMaxScore = iScore;
                pOptimalStart = pLastStart;
                pOptimalEnd = pSeed;
            } // if
        } // for

        /*
         * we need to increment both iterator but check if they extend past the end of pIn
         * (check twice because incrementing an iterator pointing to end will
         * result in undefined behaviour)
         *
         * We then move all known suboptimal seeds to pRet and swap the shared pointers.
         *
         * @note: order is significant!
         * """
         * --- Iterator validity ---
         * Iterators, pointers and references pointing to position (or first) and beyond are
         * invalidated, with all iterators, pointers and references to elements
         * before position (or first) are guaranteed to keep referring to the same
         * elements they were referring to before the call.
         * """
         */
#define KEEP_SMALLER_PART ( 0 )
#if KEEP_SMALLER_PART == 1
        if( pOptimalStart != pIn->end( ) )
            for( auto pCur = pIn->begin( ); pCur != pOptimalStart; pCur++ )
                pRet->push_back( *pCur );
#endif
        if( pOptimalEnd != pIn->end( ) )
            if( ++pOptimalEnd != pIn->end( ) )
            {
#if KEEP_SMALLER_PART == 1
                for( auto pCur = pOptimalEnd; pCur != pIn->end( ); pCur++ )
                    pRet->push_back( *pCur );
#endif
                pIn->erase( pOptimalEnd, pIn->end( ) );
            }
        if( pOptimalStart != pIn->end( ) )
            pIn->erase( pIn->begin( ), pOptimalStart );
    } // if
    // swap the shared pointers
    pRet.swap( pIn );
    /*
     * Artifact filter:
     */
    if( pRet->size( ) > 2 )
    {
        size_t uiPrePos = 0;
        size_t uiCenterPos = 1;
        while( uiCenterPos < pRet->size( ) - 1 )
        {
            Seed& rPre = ( *pRet )[ uiPrePos ];
            Seed& rCenter = ( *pRet )[ uiCenterPos ];
            Seed& rPost = ( *pRet )[ uiCenterPos + 1 ];

            int64_t iDeltaPre = rPre.start_ref( ) - (int64_t)rPre.start( );
            int64_t iDeltaCenter = rCenter.start_ref( ) - (int64_t)rCenter.start( );
            int64_t iDeltaPost = rPost.start_ref( ) - (int64_t)rPost.start( );

            int64_t iDeltaDistToPre = std::abs( iDeltaPre - iDeltaCenter );
            int64_t iDeltaDistToPost = std::abs( iDeltaPost - iDeltaCenter );

            double dDeltaDistDiff =
                std::abs( iDeltaDistToPre - iDeltaDistToPost ) * 2 / ( (double)iDeltaDistToPre + iDeltaDistToPost );

            int64_t iScoreSeed = rCenter.size() * pGlobalParams->iMatch->get( );
            int64_t iIndelPenalty = std::min((int64_t)pGlobalParams->uiSVPenalty->get( ) * 2, 
                                    (int64_t)(pGlobalParams->iGap->get( ) * 2 + 
                                        pGlobalParams->iExtend->get( ) * (iDeltaDistToPre + iDeltaDistToPost))
                                );

            if( dDeltaDistDiff < dMaxDeltaDist && (uint64_t)iDeltaDistToPre > uiMinDeltaDist && 
                iScoreSeed < iIndelPenalty * 2 )
            {
                // the filter triggered ->
                // we flag the seed to be ignored by setting it's size to zero
                rCenter.size( 0 );
                uiCenterPos++;
            } // if
            else
            {
                // the filter did not trigger move to the next triple
                uiCenterPos++;
                uiPrePos = uiCenterPos - 1;
            } // else
        } // while
    } // if

    return pRet;
} // method

inline nucSeqIndex difference( nucSeqIndex a, nucSeqIndex b )
{
    if( a > b )
        return a - b;
    return b - a;
} // function

std::shared_ptr<std::vector<std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>>>
Harmonization::linesweep( std::shared_ptr<std::vector<std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>>> pShadows,
                          const int64_t uiRStart, const double fAngle )
{
    // sort shadows (increasingly) by start coordinate of the match
    std::sort( pShadows->begin( ), pShadows->end( ),
               []( std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex> xA,
                   std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>
                       xB ) {
                   /*
                    * sort by the interval starts
                    * if two intervals start at the same point the larger one shall be treated first
                    */
                   if( std::get<1>( xA ) == std::get<1>( xB ) )
                       return ( std::get<2>( xA ) > std::get<2>( xB ) );
                   return ( std::get<1>( xA ) < std::get<1>( xB ) );
               } // lambda
    ); // sort function call

    // for(std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>& xTup : *pShadows)
    //{
    //    std::cout << std::get<0>(xTup)->start() << ", " << std::get<0>(xTup)->start_ref() << ", "
    //    << std::get<0>(xTup)->size() << std::endl;
    //}

    auto pItervalEnds = std::make_shared<std::vector<std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>>>( );
    pItervalEnds->reserve( pShadows->size( ) );

    nucSeqIndex x = 0;
    // this is the line sweeping part
    for( auto& xTup : *pShadows )
    {
        if( x < std::get<2>( xTup ) )
        {
            pItervalEnds->push_back( xTup );
            x = std::get<2>( xTup );
        } // if
        else
        {
            assert( !pItervalEnds->empty( ) );
            double fDistance = deltaDistance( *std::get<0>( xTup ), fAngle, uiRStart );
            ;
            // std::cout << "D-";
            nucSeqIndex uiPos = pItervalEnds->size( ); // uiPos is unsigned!!!
            bool bThisIsCloserToDiagonal = true;
            while( uiPos > 0 && std::get<2>( ( *pItervalEnds )[ uiPos - 1 ] ) >= std::get<2>( xTup ) )
            {
                double fDistanceOther =
                    deltaDistance( *std::get<0>( ( *pItervalEnds )[ uiPos - 1 ] ), fAngle, uiRStart );
                if( fDistanceOther <= fDistance )
                {
                    bThisIsCloserToDiagonal = false;
                    break;
                } // if
                --uiPos;
            } // while
            if( bThisIsCloserToDiagonal )
            {
                while( !pItervalEnds->empty( ) && std::get<2>( pItervalEnds->back( ) ) >= std::get<2>( xTup ) )
                    pItervalEnds->pop_back( );
                pItervalEnds->push_back( xTup );
            } // if
            // else do nothing
        } // else
    } // for
    // std::cout << pItervalEnds->size() << std::endl;
    return pItervalEnds;
} // function

std::shared_ptr<Seeds> Harmonization::harmonizeOne( std::shared_ptr<Seeds>& pSeedsIn, std::shared_ptr<NucSeq> pQuery,
                                                    std::shared_ptr<SoCPriorityQueue> pSoCIn )
{
    auto pSeeds = std::make_shared<Seeds>( );
    pSeeds->xStats = pSeedsIn->xStats;
    if( pSeedsIn->size( ) > 1 )
    {
        
        DEBUG( pSoCIn->vSoCs.push_back( pSeedsIn ); ) // DEBUG
#if USE_RANSAC == 1 // switch between ransac line angle + intercept estimation & 45deg median line
        std::vector<double> vX, vY;
        vX.reserve( pSeedsIn->size( ) );
        vY.reserve( pSeedsIn->size( ) );
        for( const auto& rSeed : *pSeedsIn )
        {
            // middlepoint of seed in plane
            vX.push_back( (double)rSeed.start_ref( ) + rSeed.size( ) / 2.0 );
            vY.push_back( (double)rSeed.start( ) + rSeed.size( ) / 2.0 );
            // start point of seed in plane
            vX.push_back( (double)rSeed.start_ref( ) );
            vY.push_back( (double)rSeed.start( ) );
            // end point of seed in plane
            vX.push_back( (double)rSeed.start_ref( ) + rSeed.size( ) );
            vY.push_back( (double)rSeed.start( ) + rSeed.size( ) );
        } // for
        /* The Mean Absolute Deviation (MAD) is later required for the threshold t */
        double fMAD = medianAbsoluteDeviation<double>( vY );
        auto xSlopeIntercept = run_ransac( vX, vY, /*pSoCIn->vIngroup.back(),*/ fMAD );

        /*
         * remove outliers
         */
        pSeedsIn->erase( std::remove_if( pSeedsIn->begin( ), pSeedsIn->end( ),
                                         [ & ]( const Seed& rS ) {
                                             return deltaDistance( rS, xSlopeIntercept.first,
                                                                   (int64_t)xSlopeIntercept.second ) > fMAD && 
                                                                   rS.size() < 50;
                                         } ),
                         pSeedsIn->end( ) );

#else
        auto rMedianSeed = ( *pSeedsIn )[ pSeedsIn->size( ) / 2 ];
        auto xSlopeIntercept = std::make_pair( 0.785398, // forty five degrees
                                               (double)rMedianSeed.start_ref( ) - (double)rMedianSeed.start( ) );
#endif

        DEBUG( pSoCIn->vSlopes.push_back( std::tan( xSlopeIntercept.first ) );
               pSoCIn->vIntercepts.push_back( xSlopeIntercept.second ); ) // DEBUG

        auto pShadows = std::make_shared<std::vector<std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>>>( );
        pShadows->reserve( pSeedsIn->size( ) );

        // get the left shadows
        for( Seeds::iterator pSeed = pSeedsIn->begin( ); pSeed != pSeedsIn->end( ); pSeed++ )
        {
            pShadows->push_back( std::make_tuple( pSeed, pSeed->start( ), pSeed->end_ref( ) ) );
            //std::cout << "1 (" << pSeed->start() << "," << pSeed->start_ref() << "," <<
            //pSeed->size() << ")" << std::endl;
        } // for

        // perform the line sweep algorithm on the left shadows
        auto pShadows2 = linesweep( pShadows, (int64_t)xSlopeIntercept.second, xSlopeIntercept.first );
        pShadows->clear( );
        pShadows->reserve( pShadows2->size( ) );

        // get the right shadows
        for( auto& xT : *pShadows2 )
        {
            pShadows->push_back(
                std::make_tuple( std::get<0>( xT ), std::get<0>( xT )->start_ref( ), std::get<0>( xT )->end( ) ) );
            //std::cout << "2 (" << std::get<0>( xT )->start_ref() << "," << std::get<0>( xT )->end() << "," <<
            //std::get<0>( xT )->size() << ")" << std::endl;
        }

        // perform the line sweep algorithm on the right shadows
        pShadows = linesweep( pShadows, (int64_t)xSlopeIntercept.second, xSlopeIntercept.first );

        pSeeds->reserve( pShadows->size( ) );

        for( auto& xT : *pShadows )
        {
            pSeeds->push_back( *std::get<0>( xT ) );
            // std::cout << "3 (" << std::get<0>( xT )->start() << "," << std::get<0>( xT )->start_ref() << "," <<
            // std::get<0>( xT )->size() << ")" << std::endl;
        }

        pSeeds->bConsistent = true;

        DEBUG( pSoCIn->vHarmSoCs.push_back( pSeeds ); ) // DEBUG

        // seeds need to be sorted for the following steps
        std::sort( pSeeds->begin( ), pSeeds->end( ),
                   []( const Seed& xA, const Seed& xB ) {
                       if( xA.start_ref( ) == xB.start_ref( ) )
                           return xA.start( ) < xB.start( );
                       return xA.start_ref( ) < xB.start_ref( );
                   } // lambda
        ); // sort function call

        DEBUG( if( !pSeeds->empty( ) ) pSoCIn->vExtractOrder.back( ).rStart = pSeeds->front( ).start_ref( );
               if( !pSeeds->empty( ) ) pSoCIn->vExtractOrder.back( ).rEnd = pSeeds->back( ).end_ref( ); ) // DEBUG

        /*
         * sometimes we have one or less seeds remaining after the cupling:
         * in these cases we simply return the center seed (judging by the delta from SoCs)
         * This increases accuracy since the aligner is guided towards the correct
         * position rather than a random seed or no seed...
         *
         * Note: returning the longest or least ambiguous seed does not make sense here;
         * In most cases where this condition triggeres we have seeds that are roughly the
         * same length and ambiguity... (e.g. 3 seeds of length 17 and two of length 16)
         */
        if( pSeeds->size( ) <= 1 )
        {
            pSeeds->clear( );
            pSeeds->push_back( ( *pSeedsIn )[ pSeedsIn->size( ) / 2 ] );
        } // if
        assert( !pSeeds->empty( ) );

    } // if
    else if( !pSeedsIn->empty( ) ) // pSeedsIn contains merely one seed
    {
        pSeeds->push_back( pSeedsIn->front( ) );
        DEBUG( pSoCIn->vExtractOrder.push_back( SoCPriorityQueue::blub( ) );
               pSoCIn->vExtractOrder.back( ).first = pSeedsIn->front( ).size( );
               pSoCIn->vExtractOrder.back( ).rStartSoC = pSeedsIn->front( ).start_ref( );
               pSoCIn->vExtractOrder.back( ).rEndSoC = pSeedsIn->front( ).end_ref( );
               pSoCIn->vIngroup.push_back( std::make_shared<Seeds>( ) );
               pSoCIn->vSlopes.push_back( 0 );
               pSoCIn->vIntercepts.push_back( 0 );
               pSoCIn->vHarmSoCs.push_back( std::make_shared<Seeds>( pSeeds ) );
               pSoCIn->vSoCs.push_back( std::make_shared<Seeds>( pSeeds ) ); ) // DEBUG
    } // else
    return pSeeds;
}
std::shared_ptr<ContainerVector<std::shared_ptr<Seeds>>>
Harmonization::execute( std::shared_ptr<SoCPriorityQueue> pSoCIn, std::shared_ptr<NucSeq> pQuery,
                        std::shared_ptr<FMIndex> pFM_index )
{
#define FILTER_1 ( 0 )
#if FILTER_1
    nucSeqIndex uiAccumulativeSeedLength = 0;
#endif

    unsigned int uiNumTries = 0;
    nucSeqIndex uiLastHarmScore = 0;
    nucSeqIndex uiBestSoCScore = 0;
    unsigned int uiSoCRepeatCounter = 0;

    auto pSoCs = std::make_shared<ContainerVector<std::shared_ptr<Seeds>>>( );

    while( !pSoCIn->empty( ) )
    {
        if( ++uiNumTries > uiMaxTries )
        {
            PRINT_BREAK_CRITERIA( std::cout << "break after " << uiNumTries << " tries." << std::endl; )
            break;
        }
        auto pSeedsIn = pSoCIn->pop( );

        DEBUG( pSoCIn->vExtractOrder.push_back( SoCPriorityQueue::blub( ) );
               pSoCIn->vExtractOrder.back( ).rStartSoC = pSeedsIn->front( ).start_ref( );
               pSoCIn->vExtractOrder.back( ).rEndSoC = pSeedsIn->front( ).end_ref( ); ) // DEBUG

        assert( !pSeedsIn->empty( ) );
        nucSeqIndex uiCurrSoCScore = 0;
        for( const auto& rSeed : *pSeedsIn )
        {
            uiCurrSoCScore += rSeed.size( );
            DEBUG( pSoCIn->vExtractOrder.back( ).rStartSoC =
                       std::min( pSoCIn->vExtractOrder.back( ).rStartSoC, rSeed.start_ref( ) );
                   pSoCIn->vExtractOrder.back( ).rEndSoC = std::max( pSoCIn->vExtractOrder.back( ).rEndSoC,
                                                                     rSeed.end_ref( ) ); ) // DEBUG
        } // for

        // Prof. Kutzners filter:
        // this merely checks weather we actually do have to do the harmonization at all
        if( bDoHeuristics && uiNumTries > uiMinTries )
        {
            if( pQuery->length( ) > uiSwitchQLen && uiSwitchQLen != 0 )
            {
                if( uiLastHarmScore > uiCurrSoCScore )
                {
                    PRINT_BREAK_CRITERIA( std::cout << "skip because of SoC score minimum" << std::endl; )
                    continue;
                } // if
            } // if

            if( uiBestSoCScore * fScoreTolerace > uiCurrSoCScore && fScoreTolerace > 0 )
            {
                PRINT_BREAK_CRITERIA( std::cout << "Break because of fast SoC drop"; )
                break;
            } // if
        } // if
        uiBestSoCScore = std::max( uiBestSoCScore, uiCurrSoCScore );
        DEBUG( pSoCIn->vExtractOrder.back( ).first = uiCurrSoCScore;
               pSoCIn->vIngroup.push_back( std::make_shared<Seeds>( ) ); )

        auto pUnharmSeedsRev = pSeedsIn->extractStrand( false );
        for( auto& rSeed : *pUnharmSeedsRev )
        {
            assert( rSeed.end( ) <= pQuery->length( ) );
            rSeed.uiPosOnReference = pFM_index->getRefSeqLength( ) - rSeed.uiPosOnReference - 1;
        }

        auto pSeedsForw = harmonizeOne( pSeedsIn, pQuery, pSoCIn );
        auto pSeedsRev = harmonizeOne( pUnharmSeedsRev, pQuery, pSoCIn );

        PRINT_BREAK_CRITERIA( if( pSoCIn->empty( ) ) std::cout << "exhausted all SoCs" << std::endl; )

        /*
         * end of FILTER
         */
        nucSeqIndex uiCurrHarmScore = 0;
        for( const auto& rSeed : *pSeedsForw )
            uiCurrHarmScore += rSeed.size( );
        for( const auto& rSeed : *pSeedsRev )
            uiCurrHarmScore += rSeed.size( );
        if( bDoHeuristics && uiNumTries > uiMinTries )
            if( uiCurrHarmScore < uiCurrHarmScoreMin )
            {
                PRINT_BREAK_CRITERIA( std::cout << "skip because of abs. harmonization score minimum" << std::endl; )
                continue;
            } // if
        if( bDoHeuristics )
        {
            if( uiCurrHarmScore < pQuery->length( ) * fCurrHarmScoreMinRel )
            {
                PRINT_BREAK_CRITERIA( std::cout << "skip because of rel. harmonization score minimum" << std::endl; )
                continue;
            } // if
        } // if
#if DEBUG_LEVEL >= 1
        std::vector<bool> vQCoverage( pQuery->length( ), false );
        pSoCIn->vExtractOrder.back( ).qCoverage = 0;
        for( const auto& rSeed : *pSeedsForw )
        {
            for( auto uiX = rSeed.start( ); uiX < rSeed.end( ); uiX++ )
                if( uiX >= 0 && uiX < pQuery->length( ) )
                    vQCoverage[ uiX ] = true;
        } // for
        for( const auto& rSeed : *pSeedsRev )
        {
            for( auto uiX = rSeed.start( ); uiX < rSeed.end( ); uiX++ )
                if( uiX >= 0 && uiX < pQuery->length( ) )
                    vQCoverage[ uiX ] = true;
        } // for
        pSoCIn->vExtractOrder.back( ).second = uiCurrHarmScore;

        for( bool b : vQCoverage )
            if( b )
                pSoCIn->vExtractOrder.back( ).qCoverage++;
#endif
        if( bDoHeuristics && uiNumTries > uiMinTries && pQuery->length( ) > uiSwitchQLen && uiSwitchQLen != 0 )
        {
            // Prof. Kutzners filter:
            if( uiLastHarmScore > uiCurrHarmScore )
            {
                PRINT_BREAK_CRITERIA( std::cout << "skip because of harmonization score dropoff" << std::endl; )
                continue;
            } // if
        } // if
        while( !pSeedsForw->empty( ) )
        {
            uiSoCRepeatCounter++;
            pSoCs->push_back( this->applyFilters( pSeedsForw ) );
        }
        while( !pSeedsRev->empty( ) )
        {
            uiSoCRepeatCounter++;
            pSoCs->push_back( this->applyFilters( pSeedsRev ) );
        }

        if( bDoHeuristics && uiNumTries > uiMinTries && pQuery->length( ) < uiSwitchQLen && uiSwitchQLen != 0 )
        {
            if( !( uiCurrHarmScore + ( pQuery->length( ) * fScoreDiffTolerance ) >= uiLastHarmScore &&
                   uiCurrHarmScore - ( pQuery->length( ) * fScoreDiffTolerance ) <= uiLastHarmScore ) )
                uiSoCRepeatCounter = 0;
            if( uiSoCRepeatCounter >= uiMaxEqualScoreLookahead && uiMaxEqualScoreLookahead != 0 )
            {
                PRINT_BREAK_CRITERIA( std::cout << "Break because of repeated harmonization"
                                                << " score for short queries " << uiMaxEqualScoreLookahead
                                                << std::endl; )
                break;
            } // else
        } // if
        else
            uiSoCRepeatCounter = 0;
        uiLastHarmScore = uiCurrHarmScore;

        // FILTER
#if FILTER_1
        if( bDoHeuristics && uiNumTries > uiMinTries )
        {
            nucSeqIndex uiAccLen = pSeedsForw->getScore( ) + pSeedsRev->getScore( );
            if( uiAccumulativeSeedLength > uiAccLen )
            {
                PRINT_BREAK_CRITERIA( std::cout << "skip because of accumulative seed length" << std::endl; )
                continue;
            } // if
            uiAccumulativeSeedLength = std::max( uiAccLen, uiAccumulativeSeedLength );
        } // if
#endif
        // FILTER END

    } // while

    if( bDoHeuristics )
    {
        for( unsigned int ui = 0; ui < uiSoCRepeatCounter && pSoCs->size( ) > uiMinTries; ui++ )
            pSoCs->pop_back( );
    } // if

    PRINT_BREAK_CRITERIA( std::cout << "computed " << pSoCs->size( ) << " SoCs." << std::endl; ) // DEBUG

    return pSoCs;
} // function

#ifdef WITH_PYTHON

void exportHarmonization( libMS::SubmoduleOrganizer& xOrganizer )
{

    exportModule<SeedLumping>( xOrganizer, "SeedLumping", []( auto&& x ) {
        x.def( "lump", &SeedLumping::lump ).def( "execute_py", &SeedLumping::execute_py );
    } );
    exportModule<SeedExtender>( xOrganizer, "SeedExtender",
                                []( auto&& x ) { x.def( "extend", &SeedExtender::extend ); } );
    exportModule<MaxExtendedToSMEM>( xOrganizer, "MaxExtendedToSMEM",
                                     []( auto&& x ) { x.def( "filter", &MaxExtendedToSMEM::filter ); } );
    exportModule<MaxExtendedToMaxSpanning>( xOrganizer, "MaxExtendedToMaxSpanning",
                                            []( auto&& x ) { x.def( "filter", &MaxExtendedToMaxSpanning::filter ); } );
    exportModule<MinLength, size_t>( xOrganizer, "MinLength",
                                     []( auto&& x ) { x.def( "filter", &MinLength::filter ); } );
    exportModule<FilterContigBorder>( xOrganizer, "FilterContigBorder" );
    exportModule<FilterOverlappingSeeds>( xOrganizer, "FilterOverlappingSeeds",
                                          []( auto&& x ) { x.def( "filter", &FilterOverlappingSeeds::filter ); } );
    exportModule<SortRemoveDuplicates>( xOrganizer, "SortRemoveDuplicates",
                                        []( auto&& x ) { x.def( "filter", &SortRemoveDuplicates::filter ); } );

    exportModule<FilterToUnique>( xOrganizer, "FilterToUnique", []( auto&& x ) {
        x.def_readwrite( "num_mm", &FilterToUnique::uiNumMissmatchesAllowed );
    } );

    exportModule<ParlindromeFilter>( xOrganizer, "ParlindromeFilter", []( auto&& x ) {
        x.def_readwrite( "palindromes", &ParlindromeFilter::pParlindromes )
            .def( "keep_palindromes", &ParlindromeFilter::keepParlindromes );
    } );

    exportModule<Harmonization>( xOrganizer, "Harmonization" );

    py::bind_vector_ext<ContainerVector<std::shared_ptr<Seeds>>, libMS::Container,
                        std::shared_ptr<ContainerVector<std::shared_ptr<Seeds>>>>( xOrganizer.container( ),
                                                                                   "ContainerVectorSeeds", "docstr" )
        .def( py::init<>( ) );

    // tell boost python that pointers of these classes can be converted implicitly
    py::implicitly_convertible<ContainerVector<std::shared_ptr<Seeds>>, libMS::Container>( );
} // function
#endif