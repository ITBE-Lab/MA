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
 */
class LastMatchingSeeds
{
    std::vector<size_t> vSeedSizes;
    std::vector<CyclicSetWithNElements<Seed>> vContent;

  public:
    LastMatchingSeeds( std::vector<size_t> vSeedSizes = {18, 30}, std::vector<size_t> vSetSizes = {2, 2} )
        : vSeedSizes( vSeedSizes )
    {
        assert( vSeedSizes.size( ) == vSetSizes.size( ) );
        vContent.reserve( vSeedSizes.size( ) );
        for( size_t uiI : vSetSizes )
            vContent.emplace_back( uiI );
    } // constructor

    void insert( Seed& rSeed )
    {
        for( size_t uiI = 1; uiI <= vSeedSizes.size( ); uiI++ )
            if( vSeedSizes[ vSeedSizes.size( ) - uiI ] > rSeed.size( ) )
            {
                vContent[ vContent.size( ) - uiI ].insert( rSeed );
                return;
            } // if
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
                                                                    std::shared_ptr<Pack> pRefSeq,
                                                                    std::shared_ptr<FMIndex> pFM_index )
{
    auto pRet = std::make_shared<ContainerVector<SvJump>>( );

    // sort seeds by query position:
    std::sort( pSegments->begin( ), pSegments->end( ),
               []( const Segment& rA, const Segment& rB ) { return rA.start( ) < rB.start( ); } );
    Seeds vSeeds;
    // avoid multiple allocations (we can only guess the actual number of seeds here)
    vSeeds.reserve( pSegments->size( ) * 2 );
    pSegments->emplaceAllEachSeeds( *pFM_index, 0, 100, 18, vSeeds, []( ) { return true; } );

    // walk over all seeds to compute
    LastMatchingSeeds xLastSeedsForward;
    for( Seed& rCurr : vSeeds )
        helperSvJumpsFromSeedsExecute( xLastSeedsForward, rCurr, true, pRet );
    LastMatchingSeeds xLastSeedsReverse;
    for( auto itRevIt = vSeeds.rbegin( ); itRevIt != vSeeds.rend( ); itRevIt++ )
        helperSvJumpsFromSeedsExecute( xLastSeedsForward, *itRevIt, false, pRet );

    return pRet;
} // method