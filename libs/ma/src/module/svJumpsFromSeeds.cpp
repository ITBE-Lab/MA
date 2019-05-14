/**
 * @file svJumsFromSeeds.cpp
 * @author Markus Schmidt
 */
#include "module/svJumpsFromSeeds.h"

using namespace libMA;

/**
 * @brief holds a set of elements with maximal size x. if more than x elements are inserted the oldest ones are removed
 */
template <typename TP_CONTENT> class CyclicSetWithNElements
{
    std::vector<TP_CONTENT> vContent;
    size_t uiCurrPos;

  public:
    CyclicSetWithNElements( size_t uiSize ) : vContent( uiSize ), uiCurrPos( 0 )
    {} // constructor

    void insert( TP_CONTENT& rEle )
    {
        vContent[ uiCurrPos % vContent.size( ) ] = rEle;
        uiCurrPos++;
    } // method

    template <typename TP_FUNCTOR> void forall( TP_FUNCTOR&& functor )
    {
        for( size_t uiI = 0; uiI < uiCurrPos && uiI < vContent.size( ); uiI++ )
            functor( vContent[ uiI ] );
    } // method
}; // class

/**
 * @brief holds x previous seeds above size y
 * @details
 * vSeedSizes determines the minimal sizes of seeds that shall be held, while the respective vSetSizes values
 * indicate how many seeds of that size or above shall be held.
 * @note vSeedSizes must be given in DESCENDING order
 */
class LastMatchingSeeds
{
    std::vector<size_t> vSeedSizes;
    std::vector<CyclicSetWithNElements<Seed>> vContent;

  public:
    LastMatchingSeeds( std::vector<size_t> vSeedSizes = {30, 0}, std::vector<size_t> vSetSizes = {2, 2} )
        : vSeedSizes( vSeedSizes )
    {
        assert( vSeedSizes.size( ) == vSetSizes.size( ) );
        vContent.reserve( vSeedSizes.size( ) );
        for( size_t uiI : vSetSizes )
            vContent.emplace_back( uiI );
        for( size_t uiI = 1; uiI < vSeedSizes.size( ); uiI++ )
            assert( vSeedSizes[ uiI - 1 ] > vSeedSizes[ uiI ] );
    } // constructor

    void insert( Seed& rSeed )
    {
        for( size_t uiI = 0; uiI < vSeedSizes.size( ); uiI++ )
            if( rSeed.size( ) >= vSeedSizes[ uiI ] )
            {
                vContent[ uiI ].insert( rSeed );
                break;
            } // if
    } // method

    template <typename TP_IT> void fill( TP_IT itBegin, TP_IT itEnd )
    {
        while( itBegin != itEnd )
        {
            insert( *itBegin );
            itBegin++;
        } // while
    } // method

    template <typename TP_FUNCTOR> void forall( TP_FUNCTOR&& functor )
    {
        for( CyclicSetWithNElements<Seed>& rX : vContent )
            rX.forall( functor );
    } // method
}; // class

void helperSvJumpsFromSeedsExecute( LastMatchingSeeds& rLastSeeds, Seed& rCurr, bool bJumpFromStart,
                                    std::shared_ptr<ContainerVector<SvJump>>& pRet )
{
    // @todo filter segment here
    rLastSeeds.forall( [&]( Seed& rLast ) //
                       {
                           if( SvJump::validJump( rLast, rCurr, bJumpFromStart ) )
                               pRet->emplace_back( rLast, rCurr, bJumpFromStart );
                       } // lambda
    ); // forall function call

    rLastSeeds.insert( rCurr );
} // function

std::shared_ptr<ContainerVector<SvJump>> SvJumpsFromSeeds::execute( std::shared_ptr<SegmentVector> pSegments,
                                                                    std::shared_ptr<Pack>
                                                                        pRefSeq,
                                                                    std::shared_ptr<FMIndex>
                                                                        pFM_index,
                                                                    std::shared_ptr<NucSeq>
                                                                        pQuery )
{
    auto pRet = std::make_shared<ContainerVector<SvJump>>( );

    // sort seeds by query position:
    std::sort( pSegments->begin( ), pSegments->end( ),
               []( const Segment& rA, const Segment& rB ) { return rA.start( ) < rB.start( ); } );
    Seeds vSeeds;
    // avoid multiple allocations (we can only guess the actual number of seeds here)
    vSeeds.reserve( pSegments->size( ) * 2 );
    pSegments->emplaceAllEachSeeds( *pFM_index, 0, uiMaxAmbiguitySv, uiMinSeedSizeSV, vSeeds, []( ) { return true; } );

    if( vSeeds.size( ) > 0 && bDoDummyJumps )
    {
        if( vSeeds.front( ).start( ) > uiMinDistDummy )
            pRet->emplace_back( vSeeds.front( ), pQuery->length( ), true );
        if( vSeeds.back( ).end( ) + uiMinDistDummy < pQuery->length( ) )
            pRet->emplace_back( vSeeds.back( ), pQuery->length( ), false );
    } // if

    // walk over all seeds to compute
    LastMatchingSeeds xLastSeedsForward;
    for( Seed& rCurr : vSeeds )
        helperSvJumpsFromSeedsExecute( xLastSeedsForward, rCurr, false, pRet );
    LastMatchingSeeds xLastSeedsReverse;
    for( auto itRevIt = vSeeds.rbegin( ); itRevIt != vSeeds.rend( ); itRevIt++ )
        helperSvJumpsFromSeedsExecute( xLastSeedsReverse, *itRevIt, true, pRet );

    return pRet;
} // method

std::shared_ptr<ContainerVector<SvJump>> SvJumpsFromSeedsPaired::execute( std::shared_ptr<SegmentVector> pSegmentsA,
                                                                          std::shared_ptr<SegmentVector>
                                                                              pSegmentsB,
                                                                          std::shared_ptr<Pack>
                                                                              pRefSeq,
                                                                          std::shared_ptr<FMIndex>
                                                                              pFM_index,
                                                                          std::shared_ptr<NucSeq>
                                                                              pQueryA,
                                                                          std::shared_ptr<NucSeq>
                                                                              pQueryB )
{
    auto pRet = std::make_shared<ContainerVector<SvJump>>( );

    // sort seeds by query position:
    std::sort( pSegmentsA->begin( ), pSegmentsA->end( ),
               []( const Segment& rA, const Segment& rB ) { return rA.start( ) < rB.start( ); } );
    Seeds vSeedsA;
    // avoid multiple allocations (we can only guess the actual number of seeds here)
    vSeedsA.reserve( pSegmentsA->size( ) * 2 );
    // emplace the seeds of the first query
    pSegmentsA->emplaceAllEachSeeds( *pFM_index, 0, uiMaxAmbiguitySv, uiMinSeedSizeSV, vSeedsA,
                                     []( ) { return true; } );
    // add the dummy jumps for the first query
    if( vSeedsA.size( ) > 0 && bDoDummyJumps )
    {
        if( vSeedsA.front( ).start( ) > uiMinDistDummy )
            pRet->emplace_back( vSeedsA.front( ), pQueryA->length( ), true );
        if( vSeedsA.back( ).end( ) + uiMinDistDummy < pQueryA->length( ) )
            pRet->emplace_back( vSeedsA.back( ), pQueryA->length( ), false );
    } // if


    std::sort( pSegmentsB->begin( ), pSegmentsB->end( ),
               []( const Segment& rA, const Segment& rB ) { return rA.start( ) < rB.start( ); } );
    Seeds vSeedsB;
    // avoid multiple allocations (we can only guess the actual number of seeds here)
    vSeedsB.reserve( pSegmentsB->size( ) * 2 );
    // emplace the seeds of the second query
    pSegmentsB->emplaceAllEachSeeds( *pFM_index, 0, uiMaxAmbiguitySv, uiMinSeedSizeSV, vSeedsB,
                                     []( ) { return true; } );
    // add the dummy jumps for the second query
    if( vSeedsB.size( ) > 0 && bDoDummyJumps )
    {
        if( vSeedsB.front( ).start( ) > uiMinDistDummy )
            pRet->emplace_back( vSeedsB.front( ), pQueryB->length( ), true );
        if( vSeedsB.back( ).end( ) + uiMinDistDummy < pQueryB->length( ) )
            pRet->emplace_back( vSeedsB.back( ), pQueryB->length( ), false );
    } // if

    // @note @todo from here on everything is quite inefficient

    // walk over all seeds to compute
    LastMatchingSeeds xLastSeedsForwardA;
    LastMatchingSeeds xLastSeedsForwardB;
    // fill
    xLastSeedsForwardA.fill( vSeedsB.begin( ), vSeedsB.end( ) );
    xLastSeedsForwardB.fill( vSeedsA.begin( ), vSeedsA.end( ) );

    LastMatchingSeeds xLastSeedsReverseA;
    LastMatchingSeeds xLastSeedsReverseB;
    // fill
    xLastSeedsForwardA.fill( vSeedsB.rbegin( ), vSeedsB.rend( ) );
    xLastSeedsForwardB.fill( vSeedsA.rbegin( ), vSeedsA.rend( ) );

    // now adjust seeds positions using the paired read information.
    for( Seed& rCurr : vSeedsA )
        rCurr.iStart += uiPairedDist;
    for( Seed& rCurr : vSeedsB )
        rCurr.iStart += uiPairedDist;


    for( Seed& rCurr : vSeedsA )
        helperSvJumpsFromSeedsExecute( xLastSeedsForwardA, rCurr, false, pRet );
    for( Seed& rCurr : vSeedsB )
        helperSvJumpsFromSeedsExecute( xLastSeedsForwardB, rCurr, false, pRet );


    for( auto itRevIt = vSeedsA.rbegin( ); itRevIt != vSeedsA.rend( ); itRevIt++ )
        helperSvJumpsFromSeedsExecute( xLastSeedsReverseA, *itRevIt, true, pRet );
    for( auto itRevIt = vSeedsB.rbegin( ); itRevIt != vSeedsB.rend( ); itRevIt++ )
        helperSvJumpsFromSeedsExecute( xLastSeedsReverseB, *itRevIt, true, pRet );

    return pRet;
} // method

#ifdef WITH_PYTHON
void exportSvJumpsFromSeeds( py::module& rxPyModuleId )
{
    exportModule<SvJumpsFromSeeds>( rxPyModuleId, "SvJumpsFromSeeds" );
    exportModule<SvJumpsFromSeedsPaired>( rxPyModuleId, "SvJumpsFromSeedsPaired" );
} // function
#endif