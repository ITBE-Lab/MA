/**
 * @file harmonization.cpp
 * @author Markus Schmidt
 */
#include "module/harmonization.h"
#if USE_RANSAC == 1
#include "sample_consensus/test_ransac.h"
#endif
using namespace libMA;

//#define PRINT_BREAK_CRITERIA(x) DEFINE(x)
#define PRINT_BREAK_CRITERIA( x )

using namespace libMA::defaults;
extern int libMA::defaults::iGap;
extern int libMA::defaults::iExtend;
extern int libMA::defaults::iMatch;
extern int libMA::defaults::iMissMatch;
// extern nucSeqIndex libMA::defaults::uiMaxGapArea;

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
        int64_t iScore = iMatch * pIn->front( ).size( );
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
            iScore += iMatch * pSeed->size( );
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
                        iScore += iMatch * ( pSeed->start_ref( ) - ( pSeed - 1 )->start_ref( ) );
                } // if
                else
                {
                    if( optimisticGapEstimation )
                        iScore += iMatch * uiGap;
                    uiGap = ( pSeed->start_ref( ) - ( pSeed - 1 )->start_ref( ) ) - uiGap;
                } // else
            } // if
            uiGap *= iExtend;
            if( uiGap > 0 )
                uiGap += iGap;
#if 0
            nucSeqIndex uiGapY = 
                pSeed->start() - ( (pSeed-1)->start() + (pSeed-1)->size() );
            if( pSeed->start() < ( (pSeed-1)->start() + (pSeed-1)->size() ) )
                uiGapY = 0;
            nucSeqIndex uiGapX = 
                pSeed->start_ref() - ((pSeed-1)->start_ref() + (pSeed-1)->size());
            if( pSeed->start_ref() < ( (pSeed-1)->start_ref() + (pSeed-1)->size() ) )
                uiGapX = 0;
            if( //check for the maximal allowed gap area
                    //uiMaxGapArea == 0 -> disabled
                    ( uiMaxGapArea > 0 && uiGapX > uiMaxGapArea && uiGapY != 0 )
                        ||
                    ( uiMaxGapArea > 0 && uiGapY > uiMaxGapArea && uiGapX != 0  )
                        ||
#else
            if( uiGap > uiSVPenalty )
                uiGap = uiSVPenalty;
            if( // check for the maximal allowed gap area
#endif
            // check for negative score
                    iScore < (int64_t)uiGap
                )
                    {
                        iScore = 0;
                        pLastStart = pSeed;
                    } // if
                    else iScore -= uiGap;
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

            if( dDeltaDistDiff < dMaxDeltaDist && (uint64_t)iDeltaDistToPre > uiMinDeltaDist )
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
#if 0
    /*
        * SV detection (insertion & deletion)
        */
    for(size_t uiPos = 0; uiPos < pRet->size() - 1; uiPos++)
    {
        Seed& rA = (*pRet)[uiPos];
        if(rA.size() == 0)
            continue;// rA is already filtered out
        size_t uiBPos = uiPos+1;

        while(uiBPos < pRet->size() && (*pRet)[uiBPos].size() == 0)
            uiBPos++;// rB is already filtered out
        if(uiBPos == pRet->size())
            break;// no more non filtered seeds in vector
        Seed& rB = (*pRet)[uiBPos];

        int64_t iX = std::max(
                static_cast<int64_t>(0),
                static_cast<int64_t>(rB.start_ref())-static_cast<int64_t>(rA.end_ref())
            );
        int64_t iY = std::max(
                static_cast<int64_t>(0),
                static_cast<int64_t>(rB.start()) - static_cast<int64_t>(rA.end())
            );

        if (iX == 0 && iY == 0)
            continue;

        double dXYRatio = std::min(iX, iY) / std::max(iY, iX);

        if( 
                iX > static_cast<int64_t>(uiMaxGapArea) ||
                iY > static_cast<int64_t>(uiMaxGapArea)
                    ||
                (dXYRatio < dMaxSVRatio && std::max(iX, iY) > iMinSVDistance)
            )
            pRet->xSvInfo.vSeedIndicesOfSVIndels.push_back(uiPos + 1);
    }// for
#endif
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

std::shared_ptr<ContainerVector<std::shared_ptr<Seeds>>>
Harmonization::execute( std::shared_ptr<SoCPriorityQueue> pSoCIn, std::shared_ptr<NucSeq> pQuery )
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

        PRINT_BREAK_CRITERIA( if( pSoCIn->empty( ) ) std::cout << "exhausted all SoCs" << std::endl; )

        auto pSeeds = std::make_shared<Seeds>( );
        pSeeds->xStats = pSeedsIn->xStats;
        if( pSeedsIn->size( ) > 1 )
        {
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
                                             [&]( const Seed& rS ) {
                                                 return deltaDistance( rS, xSlopeIntercept.first,
                                                                       (int64_t)xSlopeIntercept.second ) > fMAD;
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
                // std::cout << "(" << pSeed->start() << "," << pSeed->start_ref() << "," <<
                // pSeed->size() << ")," << std::endl;
            } // for

            // perform the line sweep algorithm on the left shadows
            auto pShadows2 = linesweep( pShadows, (int64_t)xSlopeIntercept.second, xSlopeIntercept.first );
            pShadows->clear( );
            pShadows->reserve( pShadows2->size( ) );

            // get the right shadows
            for( auto& xT : *pShadows2 )
                pShadows->push_back(
                    std::make_tuple( std::get<0>( xT ), std::get<0>( xT )->start_ref( ), std::get<0>( xT )->end( ) ) );

            // perform the line sweep algorithm on the right shadows
            pShadows = linesweep( pShadows, (int64_t)xSlopeIntercept.second, xSlopeIntercept.first );

            pSeeds->reserve( pShadows->size( ) );

            for( auto& xT : *pShadows )
                pSeeds->push_back( *std::get<0>( xT ) );

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

            DEBUG( if(!pSeeds->empty()) pSoCIn->vExtractOrder.back( ).rStart = pSeeds->front( ).start_ref( );
                   if(!pSeeds->empty()) pSoCIn->vExtractOrder.back( ).rEnd = pSeeds->back( ).end_ref( ); ) // DEBUG

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
        else // pSeedsIn contains merely one seed
        {
            pSeeds->push_back( pSeedsIn->front( ) );
            DEBUG( pSoCIn->vExtractOrder.push_back( SoCPriorityQueue::blub( ) );
                   pSoCIn->vExtractOrder.back( ).first = pSeedsIn->front( ).size( );
                   pSoCIn->vExtractOrder.back( ).rStartSoC = pSeedsIn->front( ).start_ref( );
                   pSoCIn->vExtractOrder.back( ).rEndSoC = pSeedsIn->front( ).end_ref( );
                   pSoCIn->vIngroup.push_back( std::make_shared<Seeds>( ) );
                   pSoCIn->vSlopes.push_back( 0 );
                   pSoCIn->vIntercepts.push_back( 0 );
                   pSoCIn->vHarmSoCs.push_back( std::make_shared<Seeds>( pSeeds ) ); ) // DEBUG
        } // else

        /*
         * end of FILTER
         */
        nucSeqIndex uiCurrHarmScore = 0;
        for( const auto& rSeed : *pSeeds )
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
        for( const auto& rSeed : *pSeeds )
        {
            for( auto uiX = rSeed.start( ); uiX < rSeed.end( ); uiX++ )
                vQCoverage[ uiX ] = true;
        } // for
        pSoCIn->vExtractOrder.back( ).second = uiCurrHarmScore;

        for( bool b : vQCoverage )
            if( b )
                pSoCIn->vExtractOrder.back( ).qCoverage++;
#endif
        if( bDoHeuristics && uiNumTries > uiMinTries )
        {
            //@todo uiSwitchQLen != 0 should be replaced with switch
            if( pQuery->length( ) > uiSwitchQLen && uiSwitchQLen != 0 )
            {
                // Prof. Kutzners filter:
                if( uiLastHarmScore > uiCurrHarmScore )
                {
                    PRINT_BREAK_CRITERIA( std::cout << "skip because of harmonization score dropoff" << std::endl; )
                    continue;
                } // if
            } // if
            else
            {
                if( !( uiCurrHarmScore + ( pQuery->length( ) * fScoreDiffTolerance ) >= uiLastHarmScore &&
                       uiCurrHarmScore - ( pQuery->length( ) * fScoreDiffTolerance ) <= uiLastHarmScore ) )
                    uiSoCRepeatCounter = 0;
                else if( ++uiSoCRepeatCounter >= uiMaxEqualScoreLookahead && uiMaxEqualScoreLookahead != 0 )
                {
                    // cause we haven't actually pushed the current soc yet...
                    uiSoCRepeatCounter -= 1;
                    PRINT_BREAK_CRITERIA( std::cout << "Break because of repeated harmonization"
                                                    << " score for short queries " << uiMaxEqualScoreLookahead
                                                    << std::endl; )
                    break;
                } // else
            } // else
        } // if
        uiLastHarmScore = uiCurrHarmScore;

        while( !pSeeds->empty( ) )
            pSoCs->push_back( this->applyFilters( pSeeds ) );

            // FILTER
#if FILTER_1
        if( bDoHeuristics && uiNumTries > uiMinTries )
        {
            nucSeqIndex uiAccLen = pSeeds->getScore( );
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

#ifdef BOOST_PYTHON
void exportHarmonization( )
{
    exportModule<Harmonization>( "Harmonization", []( auto&& x ) {
        x.def_readwrite( "optimistic_gap_estimation", &Harmonization::optimisticGapEstimation )
            .def_readwrite( "min_coverage", &Harmonization::fMinimalQueryCoverage )
            .def_readwrite( "tolerance", &Harmonization::fScoreTolerace )
            .def_readwrite( "max_tries", &Harmonization::uiMaxTries )
            .def_readwrite( "equal_score_lookahead", &Harmonization::uiMaxEqualScoreLookahead )
            .def_readwrite( "diff_tolerance", &Harmonization::fScoreDiffTolerance )
            .def_readwrite( "switch_q_len", &Harmonization::uiSwitchQLen )
            .def_readwrite( "do_heuristics", &Harmonization::bDoHeuristics );
    } );
} // function
#else
void exportHarmonization( py::module& rxPyModuleId )
{
    exportModule<Harmonization>( rxPyModuleId, "Harmonization", []( auto&& x ) {
        x.def_readwrite( "optimistic_gap_estimation", &Harmonization::optimisticGapEstimation )
            .def_readwrite( "min_coverage", &Harmonization::fMinimalQueryCoverage )
            .def_readwrite( "tolerance", &Harmonization::fScoreTolerace )
            .def_readwrite( "equal_score_lookahead", &Harmonization::uiMaxEqualScoreLookahead )
            .def_readwrite( "diff_tolerance", &Harmonization::fScoreDiffTolerance )
            .def_readwrite( "switch_q_len", &Harmonization::uiSwitchQLen )
            .def_readwrite( "do_heuristics", &Harmonization::bDoHeuristics );
    } );
} // function
#endif
#endif