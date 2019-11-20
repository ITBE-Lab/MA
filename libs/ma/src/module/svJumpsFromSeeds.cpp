/**
 * @file svJumsFromSeeds.cpp
 * @author Markus Schmidt
 */
#include "module/svJumpsFromSeeds.h"
#include <cmath>

using namespace libMA;

/**
 * @brief holds a set of elements with maximal size x. if more than x elements are inserted the oldest ones are removed
 */
template <typename TP_CONTENT> class CyclicSetWithNElements
{
    std::vector<TP_CONTENT> vContent;
    size_t uiCurrPos;
    size_t uiCurrSize;

  public:
    CyclicSetWithNElements( size_t uiSize ) : vContent( uiSize ), uiCurrPos( 0 ), uiCurrSize( 0 )
    {} // constructor

    void insert( TP_CONTENT& rEle )
    {
        vContent[ uiCurrPos % vContent.size( ) ] = rEle;
        uiCurrPos++;
        if( uiCurrSize < vContent.size( ) )
            uiCurrSize++;
    } // method

    void del( )
    {
        if( uiCurrSize == 0 )
            return;
        uiCurrPos++;
        uiCurrSize--;
    } // method

    // if functor returns false iteration is canceled
    template <typename TP_FUNCTOR> void forall( TP_FUNCTOR&& functor )
    {
        for( size_t uiI = 0; uiI < uiCurrSize; uiI++ )
            if( !functor( vContent[ ( vContent.size( ) + uiCurrPos - ( 1 + uiI ) ) % vContent.size( ) ] ) )
                return;
    } // method

    inline size_t size( ) const
    {
        return vContent.size( );
    } // method
}; // class

/**
 * @brief holds x previous seeds above size y
 * @details
 * vSeedSizes determines the minimal sizes of seeds that shall be held, while the respective vSetSizes values
 * indicate how many seeds of that size or above shall be held.
 * @note vSeedSizes must be given in ASCENDING order
 */
class LastMatchingSeeds
{
    size_t uiFac = 2;
    std::vector<size_t> vSeedSizes;
    std::vector<CyclicSetWithNElements<Seed>> vContent;

  public:
    LastMatchingSeeds( std::vector<size_t> vSeedSizes = {0}, std::vector<size_t> vSetSizes = {1} )
        : vSeedSizes( vSeedSizes )
    {
        assert( vSeedSizes.size( ) == vSetSizes.size( ) );
        vContent.reserve( vSeedSizes.size( ) );
        for( size_t uiI : vSetSizes )
            vContent.emplace_back( uiI * uiFac );
        for( size_t uiI = 1; uiI < vSeedSizes.size( ); uiI++ )
            assert( vSeedSizes[ uiI - 1 ] < vSeedSizes[ uiI ] );
    } // constructor

    void insert( Seed& rSeed )
    {
        size_t uiI = 0;
        for( ; uiI < vSeedSizes.size( ) - 1 && rSeed.size( ) < vSeedSizes[ uiI + 1 ]; uiI++ )
            vContent[ uiI ].del( );

        vContent[ uiI ].insert( rSeed );
    } // method

    template <typename TP_IT> void fill( TP_IT itBegin, TP_IT itEnd )
    {
        while( itBegin != itEnd )
        {
            insert( *itBegin );
            itBegin++;
        } // while
    } // method

    /// functor shall return wether the seed was usable.
    /// iteration will at most deliver the by vSetSizes given number of seeds
    /// LastMatchingSeeds stores twice the by vSetSizes given number of seeds, so half of the seeds can be unusable
    /// before LMS runs out.
    template <typename TP_FUNCTOR> void forall( TP_FUNCTOR&& functor )
    {
        for( CyclicSetWithNElements<Seed>& rX : vContent )
        {
            size_t uiToGo = rX.size( ) / uiFac;
            rX.forall( [&]( Seed& rLast ) { //
                if( functor( rLast ) )
                    uiToGo--;
                return uiToGo > 0;
            } );
        }
    } // method
}; // class

void SvJumpsFromSeeds::reseedAndMakeEdge( Seed& rLast, Seed& rCurr, bool bJumpFromStart, std::shared_ptr<NucSeq> pQuery,
                                          std::shared_ptr<Pack> pRefSeq,
                                          std::shared_ptr<ContainerVector<SvJump>>& pRet )
{
    if( !SvJump::validJump( rLast, rCurr, bJumpFromStart ) )
        return;
    pRet->emplace_back( pSelectedSetting, rLast, rCurr, bJumpFromStart, pQuery->iId );
    /* if( pRet->back( ).size( ) < (nucSeqIndex)pSelectedSetting->xMinSizeEdge->get( ) )
        pRet->pop_back( );
    else*/
    if( pRet->back( ).size( ) < (nucSeqIndex)pSelectedSetting->xMaxSizeReseed->get( ) )
    {
        nucSeqIndex uiNumSupportingNt = pRet->back( ).uiNumSupportingNt;
        nucSeqIndex uiQFrom = pRet->back( ).uiQueryFrom;
        nucSeqIndex uiQTo = pRet->back( ).uiQueryTo;
        if( uiQFrom > uiQTo )
        {
            nucSeqIndex uiT = uiQFrom;
            uiQFrom = uiQTo;
            uiQTo = uiT;
        } // if
        nucSeqIndex uiRFrom = pRet->back( ).uiFrom;
        nucSeqIndex uiRTo = pRet->back( ).uiTo;
        if( uiRFrom > uiRTo )
        {
            nucSeqIndex uiT = uiRFrom;
            uiRFrom = uiRTo;
            uiRTo = uiT;
        } // if

        if( uiQTo - uiQFrom < xHashMapSeeder.uiSeedSize || uiRTo - uiRFrom < xHashMapSeeder.uiSeedSize )
            return;

        auto pQuery2 = std::make_shared<NucSeq>( pQuery->fromTo( uiQFrom, uiQTo ) );
        auto pRef = pRefSeq->vExtract( uiRFrom, uiRTo );
        auto pSeeds = xHashMapSeeder.execute( pQuery2, pRef );
        pRef->vReverseAll( );
        pRef->vSwitchAllBasePairsToComplement( );
        auto pSeeds2 = xHashMapSeeder.execute( pQuery2, pRef );

        // fix seed positions
        for( Seed& rSeed : *pSeeds )
        {
            rSeed.uiPosOnReference += uiRFrom;
            rSeed.iStart += uiQFrom;
            assert( rSeed.end( ) < pQuery->length( ) );
        } // for
        for( Seed& rSeed : *pSeeds2 )
        {
            rSeed.bOnForwStrand = false;
            // undo reversion of reference
            rSeed.uiPosOnReference = ( uiRTo - uiRFrom ) - rSeed.uiPosOnReference - 1;
            rSeed.uiPosOnReference += uiRFrom;
            rSeed.iStart += uiQFrom;
            assert( rSeed.end( ) < pQuery->length( ) );
        } // for
        // combine seeds
        pSeeds->append( pSeeds2 );

        if( pSeeds->size( ) == 0 )
            return;

        // turn k-mers into maximally extended seeds (into max. spanning seeds even)
        pSeeds->push_back( rLast );
        pSeeds->push_back( rCurr );
        auto pLumped = xSeedLumper.execute( pSeeds );

        // compute more precise jumps
        pRet->pop_back( );
        // sort in order
        std::sort( pLumped->begin( ), pLumped->end( ),
                   [&]( Seed& rA, Seed& rB ) { return rA.start( ) < rB.start( ); } );
        size_t uiSize = pRet->size( );

        if( bJumpFromStart )
            for( size_t uiI = 1; uiI < pLumped->size( ); uiI++ )
                if( SvJump::validJump( ( *pLumped )[ uiI ], ( *pLumped )[ uiI - 1 ], true ) )
                    pRet->emplace_back( pSelectedSetting, ( *pLumped )[ uiI ], ( *pLumped )[ uiI - 1 ], true,
                                        pQuery->iId );
        if( !bJumpFromStart )
            for( size_t uiI = 1; uiI < pLumped->size( ); uiI++ )
                if( SvJump::validJump( ( *pLumped )[ uiI - 1 ], ( *pLumped )[ uiI ], false ) )
                    pRet->emplace_back( pSelectedSetting, ( *pLumped )[ uiI - 1 ], ( *pLumped )[ uiI ], false,
                                        pQuery->iId );

        // increase the support (score) for the new jumps
        for( ; uiSize < pRet->size( ); uiSize++ )
            ( *pRet )[ uiSize ].uiNumSupportingNt = uiNumSupportingNt;
    } // if
} // method


void SvJumpsFromSeeds::pickSeedsForEdge( Seeds& rSeeds, std::shared_ptr<NucSeq> pQuery, std::shared_ptr<Pack> pRefSeq,
                                         std::shared_ptr<ContainerVector<SvJump>>& pRet )
{
    for( size_t uiI = 0; uiI < rSeeds.size( ); uiI++ )
    {
        if( uiI < rSeeds.size( ) - 1 )
        {
            size_t uiBest = uiI + 1;
            double dMaxScore = 0;
            for( size_t uiJ = uiBest; uiJ < rSeeds.size( ); uiJ++ )
                if( SvJump::validJump( rSeeds[ uiI ], rSeeds[ uiJ ], false ) )
                {
                    double dDist = std::abs( ( (double)rSeeds[ uiJ ].start_ref( ) ) -
                                             ( (double)rSeeds[ uiI ].bOnForwStrand
                                                   ? rSeeds[ uiI ].end_ref( ) - 1
                                                   // @note rSeeds[ uiI ] direction is mirrored on reference if
                                                   // rSeeds[ uiI ] is on rev comp strand
                                                   : rSeeds[ uiI ].start_ref( ) - rSeeds[ uiI ].size( ) + 1 ) );
                    double dScore = rSeeds[ uiJ ].size( ) / std::log10( dDist );
                    //std::cout << uiI << " -> " << uiJ << ": " << dScore << " (" << dDist << ")" << std::endl;
                    if( dScore > dMaxScore )
                    {
                        dMaxScore = dScore;
                        uiBest = uiJ;
                    } // if
                    // if the uiI seed un seed uiJ do not overlap stop checking for the next seeds
                    // if( ( rSeeds[ uiI ].end( ) <= rSeeds[ uiJ ].start( ) ||
                    //      rSeeds[ uiJ ].end( ) <= rSeeds[ uiI ].start( ) ) &&
                    //    uiJ + 1 < rSeeds.size( ) && rSeeds[ uiJ ].start( ) != rSeeds[ uiJ + 1 ].start( ) )
                    //    break;
                } // if

            if( SvJump::validJump( rSeeds[ uiI ], rSeeds[ uiBest ], false ) )
            {
                ////std::cout << "selected: " << uiI << " -> " << uiBest << std::endl;
                reseedAndMakeEdge( rSeeds[ uiI ], rSeeds[ uiBest ], false, pQuery, pRefSeq, pRet );
            } // if
        } // if
        if( uiI > 0 )
        {
            size_t uiBest = uiI - 1;
            double dMaxScore = 0;
            for( int64_t iJ = uiBest; iJ >= 0; iJ-- )
                if( SvJump::validJump( rSeeds[ uiI ], rSeeds[ iJ ], true ) )
                {
                    double dDist = std::abs( ( (double)rSeeds[ uiI ].start_ref( ) ) -
                                             ( (double)rSeeds[ iJ ].bOnForwStrand
                                                   ? rSeeds[ iJ ].end_ref( ) - 1
                                                   // @note rSeeds[ iJ ] direction is mirrored on reference if
                                                   // rSeeds[ iJ ] is on rev comp strand
                                                   : rSeeds[ iJ ].start_ref( ) - rSeeds[ iJ ].size( ) + 1 ) );
                    double dScore = rSeeds[ iJ ].size( ) / std::log10( dDist );
                    //std::cout << uiI << " <- " << iJ << ": " << dScore << " (" << dDist << ")" << std::endl;
                    if( dScore > dMaxScore )
                    {
                        dMaxScore = dScore;
                        uiBest = (size_t)iJ;
                    } // if
                    // if the uiI seed un seed uiJ do not overlap stop checking for the next seeds
                    // if( ( rSeeds[ uiI ].end( ) <= rSeeds[ iJ ].start( ) ||
                    //      rSeeds[ iJ ].end( ) <= rSeeds[ uiI ].start( ) ) &&
                    //    iJ > 0 && rSeeds[ iJ ].end( ) != rSeeds[ iJ - 1 ].end( ) )
                    //    break;
                } // if
            if( SvJump::validJump( rSeeds[ uiI ], rSeeds[ uiBest ], true ) )
            {
                //std::cout << "selected: " << uiI << " <- " << uiBest << std::endl;
                reseedAndMakeEdge( rSeeds[ uiI ], rSeeds[ uiBest ], true, pQuery, pRefSeq, pRet );
            } // if
        } // if
    } // for
} // method

void helperSvJumpsFromSeedsExecuteOuter( const std::shared_ptr<Presetting> pSelectedSetting, Seeds& rvSeeds,
                                         std::shared_ptr<ContainerVector<SvJump>>& pRet, std::shared_ptr<Pack> pRefSeq,
                                         std::shared_ptr<NucSeq> pQuery, HashMapSeeding& rHashMapSeeder,
                                         SeedLumping& rSeedLumper, AlignedMemoryManager& rMemoryManager,
                                         NeedlemanWunsch& rNMWModule, BinarySeeding& rBinarySeeding,
                                         bool bReseed = true, bool bFromStart = true, bool bFromEnd = true );

void helperSvJumpsFromSeedsExecute( const std::shared_ptr<Presetting> pSelectedSetting, LastMatchingSeeds& rLastSeeds,
                                    Seed& rCurr, bool bJumpFromStart, std::shared_ptr<ContainerVector<SvJump>>& pRet,
                                    std::shared_ptr<Pack> pRefSeq, std::shared_ptr<NucSeq> pQuery,
                                    HashMapSeeding& rHashMapSeeder, SeedLumping& rSeedLumper,
                                    AlignedMemoryManager& rMemoryManager, NeedlemanWunsch& rNMWModule,
                                    BinarySeeding& rBinarySeeding, bool bReseed = true )
{
#if DEBUG_LEVEL >= 1 // validate seed correctness
    auto sQuery2 = pQuery->fromTo( rCurr.start( ), rCurr.end( ) );
    if( !pRefSeq->bridgingSubsection( rCurr.start_ref( ), rCurr.size( ) ) )
    {
        auto sRef2 = pRefSeq->vExtract( rCurr.start_ref( ), rCurr.end_ref( ) )->toString( );
        if( !rCurr.bOnForwStrand )
        {
            auto pRef2 = pRefSeq->vExtract( rCurr.start_ref( ) - rCurr.size( ) + 1, rCurr.start_ref( ) + 1 );
            pRef2->vReverseAll( );
            pRef2->vSwitchAllBasePairsToComplement( );
            sRef2 = pRef2->toString( );
        } // if
        if( sQuery2 != sRef2 )
            std::cerr << "(helperSvJumpsFromSeedsExecute) CRITICAL: " << sQuery2 << " != " << sRef2 << std::endl;
        assert( sQuery2 == sRef2 );
    } // if
#endif
    rLastSeeds.forall(
        [&]( Seed& rLast ) //
        {
            if( !SvJump::validJump( rLast, rCurr, bJumpFromStart ) )
                return false;
            pRet->emplace_back( pSelectedSetting, rLast, rCurr, bJumpFromStart, pQuery->iId );
            /* if( pRet->back( ).size( ) < (nucSeqIndex)pSelectedSetting->xMinSizeEdge->get( ) )
                pRet->pop_back( );
            else*/
            if( pRet->back( ).size( ) < (nucSeqIndex)pSelectedSetting->xMaxSizeReseed->get( ) && bReseed )
            {
                nucSeqIndex uiNumSupportingNt = pRet->back( ).uiNumSupportingNt;
                nucSeqIndex uiQFrom = pRet->back( ).uiQueryFrom;
                nucSeqIndex uiQTo = pRet->back( ).uiQueryTo;
                if( uiQFrom > uiQTo )
                {
                    nucSeqIndex uiT = uiQFrom;
                    uiQFrom = uiQTo;
                    uiQTo = uiT;
                } // if
                nucSeqIndex uiRFrom = pRet->back( ).uiFrom;
                nucSeqIndex uiRTo = pRet->back( ).uiTo;
                if( uiRFrom > uiRTo )
                {
                    nucSeqIndex uiT = uiRFrom;
                    uiRFrom = uiRTo;
                    uiRTo = uiT;
                } // if

                if( uiQTo - uiQFrom < rHashMapSeeder.uiSeedSize || uiRTo - uiRFrom < rHashMapSeeder.uiSeedSize )
                    return true;

                auto pQuery2 = std::make_shared<NucSeq>( pQuery->fromTo( uiQFrom, uiQTo ) );
                auto pRef = pRefSeq->vExtract( uiRFrom, uiRTo );
                auto pSeeds = rHashMapSeeder.execute( pQuery2, pRef );
                pRef->vReverseAll( );
                pRef->vSwitchAllBasePairsToComplement( );
                auto pSeeds2 = rHashMapSeeder.execute( pQuery2, pRef );

                // fix seed positions
                for( Seed& rSeed : *pSeeds )
                {
                    rSeed.uiPosOnReference += uiRFrom;
                    rSeed.iStart += uiQFrom;
                    assert( rSeed.end( ) < pQuery->length( ) );
                } // for
                for( Seed& rSeed : *pSeeds2 )
                {
                    rSeed.bOnForwStrand = false;
                    // undo reversion of reference
                    rSeed.uiPosOnReference = ( uiRTo - uiRFrom ) - rSeed.uiPosOnReference - 1;
                    rSeed.uiPosOnReference += uiRFrom;
                    rSeed.iStart += uiQFrom;
                    assert( rSeed.end( ) < pQuery->length( ) );
                } // for
                // combine seeds
                pSeeds->append( pSeeds2 );

                if( pSeeds->size( ) == 0 )
                    return true;

                // turn k-mers into maximally extended seeds (into max. spanning seeds even)
                pSeeds->push_back( rLast );
                pSeeds->push_back( rCurr );
                auto pLumped = rSeedLumper.execute( pSeeds );

                // compute more precise jumps
                pRet->pop_back( );
                // sort in order
                std::sort( pLumped->begin( ), pLumped->end( ),
                           [&]( Seed& rA, Seed& rB ) { return rA.start( ) < rB.start( ); } );
                size_t uiSize = pRet->size( );
                helperSvJumpsFromSeedsExecuteOuter( pSelectedSetting, *pLumped, pRet, pRefSeq, pQuery, rHashMapSeeder,
                                                    rSeedLumper, rMemoryManager, rNMWModule, rBinarySeeding, false,
                                                    bJumpFromStart, !bJumpFromStart );
                // increase the support (score) for the new jumps
                for( ; uiSize < pRet->size( ); uiSize++ )
                    ( *pRet )[ uiSize ].uiNumSupportingNt = uiNumSupportingNt;
            } // if
            return true;
        } // lambda
    ); // forall function call

    rLastSeeds.insert( rCurr );
} // function

void helperSvJumpsFromSeedsExecuteOuter( const std::shared_ptr<Presetting> pSelectedSetting, Seeds& rvSeeds,
                                         std::shared_ptr<ContainerVector<SvJump>>& pRet, std::shared_ptr<Pack> pRefSeq,
                                         std::shared_ptr<NucSeq> pQuery, HashMapSeeding& rHashMapSeeder,
                                         SeedLumping& rSeedLumper, AlignedMemoryManager& rMemoryManager,
                                         NeedlemanWunsch& rNMWModule, BinarySeeding& rBinarySeeding, bool bReseed,
                                         bool bFromStart, bool bFromEnd )
{
    // walk over all seeds to compute
    if( bFromEnd )
    {
        LastMatchingSeeds xLastSeedsForward;
        for( Seed& rCurr : rvSeeds )
            helperSvJumpsFromSeedsExecute( pSelectedSetting, xLastSeedsForward, rCurr, false, pRet, pRefSeq, pQuery,
                                           rHashMapSeeder, rSeedLumper, rMemoryManager, rNMWModule, rBinarySeeding,
                                           bReseed );
    } // if
    if( bFromStart )
    {
        LastMatchingSeeds xLastSeedsReverse;
        for( auto itRevIt = rvSeeds.rbegin( ); itRevIt != rvSeeds.rend( ); itRevIt++ )
            helperSvJumpsFromSeedsExecute( pSelectedSetting, xLastSeedsReverse, *itRevIt, true, pRet, pRefSeq, pQuery,
                                           rHashMapSeeder, rSeedLumper, rMemoryManager, rNMWModule, rBinarySeeding,
                                           bReseed );
    } // if
} // function

void SvJumpsFromSeeds::makeJumpsByReseedingRecursive( Seed& rLast, Seed& rNext, std::shared_ptr<NucSeq> pQuery,
                            std::shared_ptr<Pack> pRefSeq, std::shared_ptr<ContainerVector<SvJump>>& pRet )
{
    //returns a (reference pos, query pos, width [on reference], height [on query])
    auto xTuple = getPositionsForSeeds(rLast, rNext, pQuery->length());
    nucSeqIndex uiR = std::get<0>(xTuple);
    nucSeqIndex uiQ = std::get<1>(xTuple);
    nucSeqIndex uiW = std::get<2>(xTuple);
    nucSeqIndex uiH = std::get<3>(xTuple);
    if()
} // function

std::shared_ptr<ContainerVector<SvJump>> SvJumpsFromSeeds::execute( std::shared_ptr<SegmentVector> pSegments,
                                                                    std::shared_ptr<Pack>
                                                                        pRefSeq,
                                                                    std::shared_ptr<FMIndex>
                                                                        pFM_index,
                                                                    std::shared_ptr<NucSeq>
                                                                        pQuery )
{
    AlignedMemoryManager xMemoryManager;
    auto pRet = std::make_shared<ContainerVector<SvJump>>( );

    // sort seeds by query position:
    std::sort( pSegments->begin( ), pSegments->end( ),
               []( const Segment& rA, const Segment& rB ) { return rA.start( ) < rB.start( ); } );
    Seeds vSeeds;
    // avoid multiple allocations (we can only guess the actual number of seeds here)
    vSeeds.reserve( pSegments->size( ) * 2 );

    // filter ambiguous segments -> 1; don't -> 0
#if 1
    /**
     * This filters segments that occur multiple times on the reference:
     * For each of those segments extract all seeds.
     * Then pick the one seed that is closest, on the reference, to:
     *  - either the last (on query) unique seed
     *  - or the next (on query) unique seed
     * This drastically reduces the number of seeds.
     * Further, it shall help dealing with ambiguous regions by eliminating all but the seeds most likely to fit
     * into a chain of seeds
     */
    std::vector<Segment*> vTemp;
    size_t uiNumSeedsTotal = 0;
    int64_t iLastUniqueRefPos = -1;
    for( Segment& rSegment : *pSegments )
    {
        if( rSegment.size( ) < uiMinSeedSizeSV )
            continue;
        uiNumSeedsTotal += rSegment.saInterval( ).size( );

        if( rSegment.saInterval( ).size( ) == 1 )
            rSegment.forEachSeed( *pFM_index, [&]( Seed& rS ) {
                // deal with vTemp
                int64_t iCurrUniqueRefPos = (int64_t)rS.start_ref( );
                for( Segment* pSeg : vTemp )
                {
                    Seed xBest;
                    int64_t iMinDist = -1;
                    pSeg->forEachSeed( *pFM_index, [&]( Seed& rS ) {
                        int64_t iStart = (int64_t)rS.start_ref( );
                        int64_t iDist = std::abs( iStart - iCurrUniqueRefPos );
                        // if the segment is not the first segment this triggers
                        if( iLastUniqueRefPos != -1 )
                            iDist = std::min( iDist, std::abs( iStart - iLastUniqueRefPos ) );
                        if( iMinDist == -1 || iDist < iMinDist )
                        {
                            xBest = rS;
                            iMinDist = iDist;
                        } // if
                        return true;
                    } ); // forEachSeed function call
                    assert( iMinDist != -1 );
                    vSeeds.push_back( xBest );
                } // for
                vTemp.clear( );

                iLastUniqueRefPos = (int64_t)rS.end_ref( );
                vSeeds.push_back( rS );
                return true;
            } ); // forEachSeed function call \ if
        else
            vTemp.push_back( &rSegment );
    } // for
    // seeds at the end
    for( Segment* pSeg : vTemp )
    {
        Seed xBest;
        int64_t iMinDist = -1;
        pSeg->forEachSeed( *pFM_index, [&]( Seed& rS ) {
            int64_t iStart = (int64_t)rS.start_ref( );
            int64_t iDist = std::abs( iStart - iLastUniqueRefPos );
            if( iMinDist == -1 || iDist < iMinDist )
            {
                xBest = rS;
                iMinDist = iDist;
            } // if
            return true;
        } ); // forEachSeed function call
        assert( iMinDist != -1 );
        vSeeds.push_back( xBest );
    } // for

    {
        std::lock_guard<std::mutex> xGuard( xLock );
        uiNumSeedsEliminatedAmbiguityFilter += uiNumSeedsTotal - vSeeds.size( );
        uiNumSeedsKeptAmbiguityFilter += vSeeds.size( );
    } // scope xGuard

#else
    pSegments->emplaceAllEachSeeds( *pFM_index, pQuery->length( ), uiMaxAmbiguitySv, uiMinSeedSizeSV, vSeeds,
                                    [&]( ) { return true; } );
#endif

    // this auto locks using the DB; seeds need to be sorted by query position
    xCoverageInserter.insert( vSeeds, pQuery->length( ) );

    if( vSeeds.size( ) > 0 && bDoDummyJumps )
    {
        if( vSeeds.front( ).start( ) > uiMinDistDummy )
            pRet->emplace_back( pSelectedSetting, vSeeds.front( ), pQuery->length( ), true, pQuery->iId );
        if( vSeeds.back( ).end( ) + uiMinDistDummy < pQuery->length( ) )
            pRet->emplace_back( pSelectedSetting, vSeeds.back( ), pQuery->length( ), false, pQuery->iId );
    } // if

#if DEBUG_LEVEL > 0
    //std::cout << "Seeds:" << std::endl;
    //size_t uiI = 0;
    //for( auto xSeed : vSeeds )
    //{
    //    xSeed.uiId = uiI;
    //    std::cout << uiI++ << ": " << xSeed.start( ) << ", " << xSeed.start_ref( ) << ", " << xSeed.size( ) << ", "
    //              << ( xSeed.bOnForwStrand ? "frow" : "rev" ) << std::endl;
    //} // for
#endif

    helperSvJumpsFromSeedsExecuteOuter( pSelectedSetting, vSeeds, pRet, pRefSeq, pQuery, xHashMapSeeder, xSeedLumper,
                                        xMemoryManager, xNMWModule, xBinarySeeding );
    //std::cout << "Edges:" << std::endl;
    //pickSeedsForEdge( vSeeds, pQuery, pRefSeq, pRet );

    return pRet;
} // method


#ifdef WITH_PYTHON
void exportSvJumpsFromSeeds( py::module& rxPyModuleId )
{
    exportModule<SvJumpsFromSeeds, int64_t, std::shared_ptr<SV_DB>, std::shared_ptr<Pack>>(
        rxPyModuleId, "SvJumpsFromSeeds", []( auto&& x ) { x.def( "commit", &SvJumpsFromSeeds::commit ); } );
} // function
#endif