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
    LastMatchingSeeds( std::vector<size_t> vSeedSizes = {0, 30}, std::vector<size_t> vSetSizes = {1, 1} )
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

void helperSvJumpsFromSeedsExecuteOuter( const std::shared_ptr<Presetting> pSelectedSetting,
                                         Seeds& rvSeeds, std::shared_ptr<ContainerVector<SvJump>>& pRet,
                                         std::shared_ptr<Pack> pRefSeq, std::shared_ptr<NucSeq> pQuery,
                                         HashMapSeeding& rHashMapSeeder, SeedLumping& rSeedLumper, bool bReseed=true );

void helperSvJumpsFromSeedsExecute( const std::shared_ptr<Presetting> pSelectedSetting, LastMatchingSeeds& rLastSeeds,
                                    Seed& rCurr, bool bJumpFromStart, std::shared_ptr<ContainerVector<SvJump>>& pRet,
                                    std::shared_ptr<Pack> pRefSeq, std::shared_ptr<NucSeq> pQuery,
                                    HashMapSeeding& rHashMapSeeder, SeedLumping& rSeedLumper, bool bReseed=true )
{
    rLastSeeds.forall(
        [&]( Seed& rLast ) //
        {
            if( SvJump::validJump( rLast, rCurr, bJumpFromStart ) )
            {
                pRet->emplace_back( pSelectedSetting, rLast, rCurr, bJumpFromStart );
                if( pRet->back( ).size( ) < (nucSeqIndex)pSelectedSetting->xMaxSizeReseed->get( ) &&
                    rLast.bOnForwStrand == rCurr.bOnForwStrand && bReseed ) // @todo make this work on inversions...
                {
                    // trigger reseeding @todo
                    nucSeqIndex uiQFrom = pRet->back( ).uiQueryFrom;
                    nucSeqIndex uiQTo = pRet->back( ).uiQueryTo;
                    if(uiQFrom > uiQTo)
                    {
                        nucSeqIndex uiT = uiQFrom;
                        uiQFrom = uiQTo;
                        uiQTo = uiT;
                    } // if
                    nucSeqIndex uiRFrom = pRet->back( ).uiFrom;
                    nucSeqIndex uiRTo = pRet->back( ).uiTo;
                    if(uiRFrom > uiRTo)
                    {
                        nucSeqIndex uiT = uiRFrom;
                        uiRFrom = uiRTo;
                        uiRTo = uiT;
                    } // if

                    auto pSeeds = rHashMapSeeder.execute(
                        std::make_shared<NucSeq>(
                            rLast.bOnForwStrand
                                ? pQuery->fromTo( uiQFrom, uiQTo )
                                : pQuery->fromToComplement( pQuery->length( ) - uiQTo, pQuery->length( ) - uiQFrom ) ),
                        pRefSeq->vExtract( uiRFrom, uiRTo ) );

                    if( pSeeds->size( ) == 0 )
                        return true;

                    for( Seed& rSeed : *pSeeds )
                    {
                        rSeed.bOnForwStrand = rLast.bOnForwStrand;
                        rSeed.uiPosOnReference += uiRFrom;
                        rSeed.iStart += uiQFrom;
                    } // for
                    // turn k-mers into maximally extended seeds
                    auto pLumped = rSeedLumper.execute( pSeeds );
                    pLumped->push_back( rLast );
                    pLumped->push_back( rCurr );

                    // compute more precise jumps
                    pRet->pop_back( ); // @note keep?
                    // sort in order
                    std::sort( pLumped->begin( ), pLumped->end( ),
                               [&]( Seed& rA, Seed& rB ) { return rA.start( ) < rB.start( ); } );
                    helperSvJumpsFromSeedsExecuteOuter( pSelectedSetting, *pLumped, pRet, pRefSeq, pQuery,
                                                        rHashMapSeeder, rSeedLumper, false );

                } // if
                return true;
            }
            return false;
        } // lambda
    ); // forall function call

    rLastSeeds.insert( rCurr );
} // function

void helperSvJumpsFromSeedsExecuteOuter( const std::shared_ptr<Presetting> pSelectedSetting,
                                         Seeds& rvSeeds, std::shared_ptr<ContainerVector<SvJump>>& pRet,
                                         std::shared_ptr<Pack> pRefSeq, std::shared_ptr<NucSeq> pQuery,
                                         HashMapSeeding& rHashMapSeeder, SeedLumping& rSeedLumper, bool bReseed )
{
    // walk over all seeds to compute
    LastMatchingSeeds xLastSeedsForward;
    for( Seed& rCurr : rvSeeds )
        helperSvJumpsFromSeedsExecute( pSelectedSetting, xLastSeedsForward, rCurr, false, pRet, pRefSeq, pQuery,
                                       rHashMapSeeder, rSeedLumper, bReseed );
    LastMatchingSeeds xLastSeedsReverse;
    for( auto itRevIt = rvSeeds.rbegin( ); itRevIt != rvSeeds.rend( ); itRevIt++ )
        helperSvJumpsFromSeedsExecute( pSelectedSetting, xLastSeedsReverse, *itRevIt, true, pRet, pRefSeq, pQuery,
                                       rHashMapSeeder, rSeedLumper, bReseed );
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
    pSegments->emplaceAllEachSeeds( *pFM_index, 0, uiMaxAmbiguitySv, uiMinSeedSizeSV, vSeeds, [&]( ) { return true; } );

    if( vSeeds.size( ) > 0 && bDoDummyJumps )
    {
        if( vSeeds.front( ).start( ) > uiMinDistDummy )
            pRet->emplace_back( pSelectedSetting, vSeeds.front( ), pQuery->length( ), true );
        if( vSeeds.back( ).end( ) + uiMinDistDummy < pQuery->length( ) )
            pRet->emplace_back( pSelectedSetting, vSeeds.back( ), pQuery->length( ), false );
    } // if

    helperSvJumpsFromSeedsExecuteOuter( pSelectedSetting, vSeeds, pRet, pRefSeq, pQuery, xHashMapSeeder, xSeedLumper );
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
            pRet->emplace_back( pSelectedSetting, vSeedsA.front( ), pQueryA->length( ), true );
        if( vSeedsA.back( ).end( ) + uiMinDistDummy < pQueryA->length( ) )
            pRet->emplace_back( pSelectedSetting, vSeedsA.back( ), pQueryA->length( ), false );
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
            pRet->emplace_back( pSelectedSetting, vSeedsB.front( ), pQueryB->length( ), true );
        if( vSeedsB.back( ).end( ) + uiMinDistDummy < pQueryB->length( ) )
            pRet->emplace_back( pSelectedSetting, vSeedsB.back( ), pQueryB->length( ), false );
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
    xLastSeedsReverseA.fill( vSeedsB.rbegin( ), vSeedsB.rend( ) );
    xLastSeedsReverseB.fill( vSeedsA.rbegin( ), vSeedsA.rend( ) );


    // now adjust seeds positions using the paired read information.
    for( Seed& rCurr : vSeedsA )
        rCurr.iStart += uiPairedDist;
    for( Seed& rCurr : vSeedsB )
        rCurr.iStart += uiPairedDist;


    for( Seed& rCurr : vSeedsA )
        helperSvJumpsFromSeedsExecute( pSelectedSetting, xLastSeedsForwardA, rCurr, false, pRet, pRefSeq, pQueryA,
                                       xHashMapSeeder, xSeedLumper );
    for( Seed& rCurr : vSeedsB )
        helperSvJumpsFromSeedsExecute( pSelectedSetting, xLastSeedsForwardB, rCurr, false, pRet, pRefSeq, pQueryB,
                                       xHashMapSeeder, xSeedLumper );


    for( auto itRevIt = vSeedsA.rbegin( ); itRevIt != vSeedsA.rend( ); itRevIt++ )
        helperSvJumpsFromSeedsExecute( pSelectedSetting, xLastSeedsReverseA, *itRevIt, true, pRet, pRefSeq, pQueryA,
                                       xHashMapSeeder, xSeedLumper );
    for( auto itRevIt = vSeedsB.rbegin( ); itRevIt != vSeedsB.rend( ); itRevIt++ )
        helperSvJumpsFromSeedsExecute( pSelectedSetting, xLastSeedsReverseB, *itRevIt, true, pRet, pRefSeq, pQueryB,
                                       xHashMapSeeder, xSeedLumper );

    return pRet;
} // method

#ifdef WITH_PYTHON
void exportSvJumpsFromSeeds( py::module& rxPyModuleId )
{
    exportModule<SvJumpsFromSeeds>( rxPyModuleId, "SvJumpsFromSeeds" );
    exportModule<SvJumpsFromSeedsPaired>( rxPyModuleId, "SvJumpsFromSeedsPaired" );
} // function
#endif