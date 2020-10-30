/**
 * @file svJumpsFromSeeds.h
 * @brief Implements a way to compute SV-jumps from seeds.
 * @author Markus Schmidt
 */
#pragma once

#include "ma/container/segment.h"
#include "ma/module/binarySeeding.h"
#include "ma/module/harmonization.h"
#include "ma/module/hashMapSeeding.h"
#include "ma/module/needlemanWunsch.h"
#include "ms/module/module.h"
#include "msv/container/svJump.h"
#include "msv/container/sv_db/tables/svJump.h"
#include "msv/module/abstractFilter.h"

namespace libMSV
{
class PerfectMatch;

#define complement( x ) ( uint8_t ) NucSeq::nucleotideComplement( x )

using RectPair = std::pair<geom::Rectangle<nucSeqIndex>, geom::Rectangle<nucSeqIndex>>;


inline bool operator<( const geom::Interval<nucSeqIndex>& xA, const geom::Interval<nucSeqIndex>& xB )
{
    if( xA.iStart != xB.iStart )
        return xA.iStart < xB.iStart;
    return xA.iSize < xB.iSize;
} // function
inline bool operator!=( const geom::Interval<nucSeqIndex>& xA, const geom::Interval<nucSeqIndex>& xB )
{
    if( xA.iStart != xB.iStart )
        return true;
    return xA.iSize != xB.iSize;
} // function

inline bool operator!=( const geom::Rectangle<nucSeqIndex>& xA, const geom::Rectangle<nucSeqIndex>& xB )
{
    if( xA.xXAxis != xB.xXAxis )
        return true;
    return xA.xYAxis != xB.xYAxis;
} // function

struct RectComp
{
    bool operator( )( const geom::Rectangle<nucSeqIndex>& xA, const geom::Rectangle<nucSeqIndex>& xB ) const
    {
        if( xA.xXAxis != xB.xXAxis )
            return xA.xXAxis < xB.xXAxis;
        return xA.xYAxis < xB.xYAxis;
    } // function
}; // struct


#define SV_JUMP_FROM_SEED_DEBUG_PRINT 0

/**
 * @brief Computes Sv-Jumps from a given seed set
 * @note WARNING: DO USE EACH INSTANCE OF THIS MODULE ONLY ONCE IN THE COMPUTATIONAL GRAPH
 */
class SvJumpsFromSeeds
    : public libMS::Module<libMS::ContainerVector<SvJump>, false, SegmentVector, Pack, FMIndex, NucSeq>
{
  public:
    const size_t uiMinSeedSizeSV;
    const size_t uiMaxAmbiguitySv;
    const int64_t iMaxSizeReseed;
    const bool bDoDummyJumps;
    const size_t uiMinDistDummy;
    const size_t uiMaxDistDummy;
    const size_t uiMinSizeJump;
    SeedLumping xSeedLumper;
    SortRemoveDuplicates xRemoveDuplicates;
    MaxExtendedToSMEM xToSMEM;
    FilterOverlappingSeeds xFilterOverlapping;
    NeedlemanWunsch xNW;
    ParlindromeFilter xParlindromeFilter;
    double dMaxSequenceSimilarity = 0.2;

    std::mutex xLock;
    size_t uiNumSeedsEliminatedAmbiguityFilter = 0;
    size_t uiNumSeedsKeptAmbiguityFilter = 0;

    /// used to indicate that there is no seed for on of the parameters in the recursive call.
    Seed xDummySeed;

    double dExtraSeedingAreaFactor = 1.5;
    // depends on sequencer technique
    // lower -> less seed noise (via larger min size for seeds); worse breakpoint recognition
    // higher -> more seed noise; better breakpoint recognition (via smaller seeds)
    double dProbabilityForRandomMatch = 0.01;

    /**
     * @brief Initialize a SvJumpsFromSeeds Module
     */
    SvJumpsFromSeeds( const ParameterSetManager& rParameters, std::shared_ptr<Pack> pRefSeq )
        : uiMinSeedSizeSV( rParameters.getSelected( )->xMinSeedSizeSV->get( ) ),
          uiMaxAmbiguitySv( rParameters.getSelected( )->xMaxAmbiguitySv->get( ) ),
          iMaxSizeReseed( (int64_t)rParameters.getSelected( )->xMaxSizeReseed->get( ) ),
          bDoDummyJumps( rParameters.getSelected( )->xDoDummyJumps->get( ) ),
          uiMinDistDummy( rParameters.getSelected( )->xMinDistDummy->get( ) ),
          uiMaxDistDummy( rParameters.getSelected( )->xMaxDistDummy->get( ) ),
          uiMinSizeJump( (size_t)rParameters.getSelected( )->xMinSizeEdge->get( ) ),
          xSeedLumper( rParameters ),
          xRemoveDuplicates( rParameters ),
          xToSMEM( rParameters ),
          xFilterOverlapping( rParameters ),
          xNW( rParameters ),
          xParlindromeFilter( rParameters )
    {} // constructor

    ~SvJumpsFromSeeds( )
    {
        size_t uiTotal = uiNumSeedsKeptAmbiguityFilter + uiNumSeedsEliminatedAmbiguityFilter;
        if( uiTotal > 0 )
            std::cout << "~SvJumpsFromSeeds: ambiguity filter kept and eliminated " << uiNumSeedsKeptAmbiguityFilter
                      << " and " << uiNumSeedsEliminatedAmbiguityFilter << " seeds respectively. " << std::endl
                      << "\tThats " << ( 1000 * uiNumSeedsKeptAmbiguityFilter / uiTotal ) / 10.0 << "% and "
                      << ( 1000 * uiNumSeedsEliminatedAmbiguityFilter / uiTotal ) / 10.0 << "% respectively."
                      << std::endl;
    } // destructor

    int64_t delta_dist( const Seed& xA, const Seed& xB )
    {
        int64_t iDeltaA = ( xA.start_ref( ) - (int64_t)xA.start( ) );
        int64_t iDeltaB = ( xB.start_ref( ) - (int64_t)xB.start( ) );

        return std::abs( iDeltaA - iDeltaB );
    }

    double overlap( const Seed& xA, const Seed& xB )
    {
        nucSeqIndex uiOverlap =
            std::max( 0, ( (int)std::min( xA.end( ), xB.end( ) ) ) - ( (int)std::max( xA.start( ), xB.start( ) ) ) );

        nucSeqIndex uiMinSize = std::min( xA.size( ), xB.size( ) );
        if( uiMinSize == 0 )
            return 1;

        return uiOverlap / (double)uiMinSize;
    }


    void forMatchingSeeds( std::shared_ptr<Seeds> pSeeds, std::function<void( Seed&, Seed& )> fOut )
    {
        for( size_t uiI = 0; uiI < pSeeds->size( ); uiI++ )
        {
            if( ( *pSeeds )[ uiI ].size( ) == 0 ) // currently at duplicate seed.
                continue;

            size_t uiJ = uiI + 1;

            // skip seeds that overlap with uiI more the 95% on query
            while( uiJ < pSeeds->size( ) &&
                   ( ( ( *pSeeds )[ uiJ ] ).size( ) == 0 || overlap( ( *pSeeds )[ uiI ], ( *pSeeds )[ uiJ ] ) > 0.95 ) )
                uiJ++;

            // uiJ now points to the first seeds that follows uiI and is not overlapping
            size_t uiK = uiJ;
            // create entry from uiI to all seeds overlapping more than 95% with uiJ (incl. itself)
            while( uiK < pSeeds->size( ) &&
                   ( ( ( *pSeeds )[ uiK ] ).size( ) == 0 || overlap( ( *pSeeds )[ uiK ], ( *pSeeds )[ uiJ ] ) > 0.95 ) )
            {
                if( ( *pSeeds )[ uiK ].size( ) != 0 ) // make sure uiK was not deleted earlier
                    fOut( ( *pSeeds )[ uiI ], ( *pSeeds )[ uiK ] );
                ++uiK;
            } // while
        } // for
    } // method

    void collectRectangles( std::shared_ptr<Seeds> pSeeds,
                            std::shared_ptr<NucSeq>
                                pQuery,
                            std::shared_ptr<Pack>
                                pRefSeq,
                            HelperRetVal* pOutExtra,
                            std::function<void( RectPair )>
                                fOut )
    {
        // first seed must be first due to the sorting
        // just skip over "deleted" (size 0) seeds here
        auto xItFirst = pSeeds->begin( );
        while( xItFirst != pSeeds->end( ) && xItFirst->size( ) == 0 )
            xItFirst++;
        if( xItFirst != pSeeds->end( ) )
            fOut( getPositionsForSeeds( xDummySeed, *xItFirst, 0, pQuery->length( ), pRefSeq ) );

        // have to look through all seeds here: with the sorting there is no guarantee where the last seed is
        // @todo this should be fixed (interval sorting?)...
        auto pLast = &pSeeds->front( );
        for( auto& rSeed : *pSeeds )
            if( rSeed.size( ) != 0 && rSeed.end( ) > pLast->end( ) )
                pLast = &rSeed;
        if( pLast->size( ) != 0 )
            fOut( getPositionsForSeeds( *pLast, xDummySeed, 0, pQuery->length( ), pRefSeq ) );

        forMatchingSeeds( pSeeds, [ & ]( Seed& rA, Seed& rB ) {
#if SV_JUMP_FROM_SEED_DEBUG_PRINT
            std::cout << "collectRectangles 1 " << rA.start( ) << ", " << rA.start_ref( ) << ", " << rA.size( )
                      << ( rA.bOnForwStrand ? " forw" : " rev" ) << std::endl;
            std::cout << "collectRectangles 2 " << rB.start( ) << ", " << rB.start_ref( ) << ", " << rB.size( )
                      << ( rB.bOnForwStrand ? " forw" : " rev" ) << std::endl;
#endif
            fOut( getPositionsForSeeds( rA, rB, 0, pQuery->length( ), pRefSeq ) );
        } );
    } // method

    void markDuplicates( std::shared_ptr<Seeds> pSeeds )
    {
        Seed* pLast = &xDummySeed;
        for( auto& xSeed : *pSeeds )
        {
            if( pLast->bOnForwStrand == xSeed.bOnForwStrand && pLast->start( ) == xSeed.start( ) &&
                pLast->size( ) == xSeed.size( ) && pLast->start_ref( ) == xSeed.start_ref( ) )
                xSeed.size( 0 ); // mark as duplicate
            else
                pLast = &xSeed;
        } // for
    } // method

    void eraseMarked( std::shared_ptr<Seeds> pSeeds )
    {
        pSeeds->erase(
            std::remove_if( pSeeds->begin( ), pSeeds->end( ), [ & ]( const Seed& xS ) { return xS.size( ) == 0; } ),
            pSeeds->end( ) );
    } // method


    /**
     * @note not threadsave if pOutExtra != nullptr
     */
    std::shared_ptr<Seeds> reseed( std::shared_ptr<Seeds> pSeeds,
                                   std::shared_ptr<NucSeq>
                                       pQuery,
                                   std::shared_ptr<Pack>
                                       pRefSeq,
                                   HelperRetVal* pOutExtra )
    {
        std::set<geom::Rectangle<nucSeqIndex>, RectComp> xReseededRectangles;
        xReseededRectangles.emplace( 0, 0, 0, 0 ); // emplace an empty rectangle so that we never reseed an empty one
        auto pRet = std::make_shared<Seeds>( );
        pRet->append( pSeeds );
        // push original seeds
        std::set<Seed, SeedCmp> xOutExtraDup;
        size_t uiOutExtraSeedsBefore = 0;
        if( pOutExtra != nullptr )
        {
            uiOutExtraSeedsBefore = pOutExtra->pSeeds->size( );
            for( auto& rSeed : *pSeeds )
            {
                if( xOutExtraDup.count( rSeed ) > 0 )
                    continue;
                xOutExtraDup.insert( rSeed );
                pOutExtra->pSeeds->push_back( rSeed );
                pOutExtra->vLayerOfSeeds.push_back( 0 );
                pOutExtra->vSocIds.push_back( pOutExtra->uiCurrSocID );
                pOutExtra->vParlindromeSeed.push_back( false );
            } // for
        } // if

        size_t uiLayer = 1;
        const size_t uiLayerLimit = 1000;
        while( uiLayer < uiLayerLimit ) // @todo emergency limit -> adjust seed sizes and rectangle size...
        {
#if SV_JUMP_FROM_SEED_DEBUG_PRINT
            std::cout << " --- Layer " << uiLayer << " --- " << std::endl;
#endif
            // sort seeds and remove duplicates (from reseeding)
            eraseMarked( pRet );
            SeedCmp xSort;
            std::sort( pRet->begin( ), pRet->end( ), xSort );
            markDuplicates( pRet );
#if SV_JUMP_FROM_SEED_DEBUG_PRINT
            for( auto& xSeed : *pRet )
                std::cout << "reseed seed " << xSeed.start( ) << ", " << xSeed.start_ref( ) << ", " << xSeed.size( )
                          << ( xSeed.bOnForwStrand ? " forw" : " rev" ) << std::endl;
#endif
            // compute new reseeding rectangles
            std::vector<geom::Rectangle<nucSeqIndex>> vNewRects;
            collectRectangles( pRet, pQuery, pRefSeq, pOutExtra, [ & ]( RectPair xNew ) {
                if( xReseededRectangles.count( xNew.first ) == 0 )
                    vNewRects.push_back( xNew.first );
                if( xReseededRectangles.count( xNew.second ) == 0 )
                    vNewRects.push_back( xNew.second );
            } );
            if( vNewRects.empty( ) )
            {
                eraseMarked( pRet );
                break;
            } // if
            // compute new seeds
            for( auto& xRect : vNewRects )
            {
#if SV_JUMP_FROM_SEED_DEBUG_PRINT
                std::cout << "reseed rect 1 " << xRect.xXAxis.start( ) << " - " << xRect.xXAxis.end( ) << ", "
                          << xRect.xYAxis.start( ) << " - " << xRect.xYAxis.end( ) << std::endl;
#endif
                xReseededRectangles.insert( xRect );
                if( pOutExtra != nullptr )
                    xParlindromeFilter.keepParlindromes( );
                auto pNewSeeds = computeSeeds( xRect, pQuery, pRefSeq, pOutExtra );
#if SV_JUMP_FROM_SEED_DEBUG_PRINT
                for( auto& xSeed : *pNewSeeds )
                    std::cout << "pNewSeeds " << xSeed.start( ) << ", " << xSeed.start_ref( ) << ", " << xSeed.size( )
                              << ( xSeed.bOnForwStrand ? " forw" : " rev" ) << std::endl;
#endif
                if( pOutExtra != nullptr )
                {
                    pOutExtra->vRectangles.push_back( xRect );
                    pOutExtra->vRectangleLayers.push_back( uiLayer );
                    pOutExtra->vRectangleFillPercentage.push_back( rectFillPercentage( pNewSeeds, xRect ) );
                } // if
                auto pParlindromeFiltered = xParlindromeFilter.execute( pNewSeeds );
#if SV_JUMP_FROM_SEED_DEBUG_PRINT
                for( auto& xSeed : *pParlindromeFiltered )
                    std::cout << "pParlindromeFiltered " << xSeed.start( ) << ", " << xSeed.start_ref( ) << ", "
                              << xSeed.size( ) << ( xSeed.bOnForwStrand ? " forw" : " rev" ) << std::endl;
#endif
                pRet->append( pParlindromeFiltered );
                if( pOutExtra != nullptr )
                {
                    for( auto& rSeed : *pParlindromeFiltered )
                    {
                        if( xOutExtraDup.count( rSeed ) > 0 )
                            continue;
                        xOutExtraDup.insert( rSeed );
                        pOutExtra->pSeeds->push_back( rSeed );
                        pOutExtra->vSocIds.push_back( pOutExtra->uiCurrSocID );
                        pOutExtra->vLayerOfSeeds.push_back( uiLayer );
                        pOutExtra->vParlindromeSeed.push_back( false );
                    } // for
                    for( auto& rSeed : *xParlindromeFilter.pParlindromes )
                    {
                        if( xOutExtraDup.count( rSeed ) > 0 )
                            continue;
                        xOutExtraDup.insert( rSeed );
                        pOutExtra->pSeeds->push_back( rSeed );
                        pOutExtra->vSocIds.push_back( pOutExtra->uiCurrSocID );
                        pOutExtra->vLayerOfSeeds.push_back( uiLayer );
                        pOutExtra->vParlindromeSeed.push_back( true );
                    } // for
                } // if
            } // for
            uiLayer++;
        } // while
        if( uiLayer == uiLayerLimit )
            std::cerr << "WARNING hit reseeding layer limit" << std::endl;

        // std::cout << uiLayer << " " << pRet->size( ) << std::endl;
        auto pFilteredRet = /*xFilterOverlapping.execute( xToSMEM.execute(*/ pRet /*) )*/;

        nucSeqIndex uiNumNtTotal = 0;
        for( auto& rSeed : *pFilteredRet )
            uiNumNtTotal += rSeed.size( );
        for( auto& rSeed : *pFilteredRet )
            rSeed.uiSoCNt = uiNumNtTotal;
        if( pOutExtra != nullptr )
        {
            // pOutExtra->computeOverlappingSeedsVector( pFilteredRet );
            for( size_t uiJ = uiOutExtraSeedsBefore; uiJ < pOutExtra->pSeeds->size( ); uiJ++ )
                ( *pOutExtra->pSeeds )[ uiJ ].uiSoCNt = uiNumNtTotal;
        } // if

        return pFilteredRet;
    } // mehtod

    std::shared_ptr<libMS::ContainerVector<SvJump>> computeJumps( std::shared_ptr<Seeds> pSeeds,
                                                                  std::shared_ptr<NucSeq>
                                                                      pQuery,
                                                                  std::shared_ptr<Pack>
                                                                      pRefSeq,
                                                                  HelperRetVal* pOutExtra )
    {
        std::sort( pSeeds->begin( ), pSeeds->end( ),
                   []( const Seed& rA, const Seed& rB ) { return rA.start( ) < rB.start( ); } );
        auto pRet = std::make_shared<libMS::ContainerVector<SvJump>>( );
        std::set<std::pair<Seed*, Seed*>> xExistingPairs;
        forMatchingSeeds( pSeeds, [ & ]( Seed& rA, Seed& rB ) {
            assert( rA.size( ) != 0 );
            assert( rB.size( ) != 0 );
            // filter out all duplicates
            if( xExistingPairs.count( std::make_pair( &rA, &rB ) ) != 0 )
                return;
            xExistingPairs.emplace( &rA, &rB );
            if( pOutExtra != nullptr )
                pOutExtra->vJumpSeeds.emplace_back( rA, rB );

            pRet->emplace_back( rA, rB, pQuery->iId );
            // remove jump again if it is too short (max of query and ref size)
            if( pRet->back( ).size( ) < uiMinSizeJump )
                pRet->pop_back( );
        } );
        // dummy jumps for first and last seed
        if( bDoDummyJumps )
        {
            auto xItFirst = pSeeds->begin( );
            while( xItFirst != pSeeds->end( ) && xItFirst->size( ) == 0 )
                xItFirst++;
            if( xItFirst != pSeeds->end( ) && xItFirst->start( ) > uiMinDistDummy )
            {
                pRet->emplace_back( *xItFirst, pQuery->length( ), true, pQuery->iId, uiMaxDistDummy );
                if( pOutExtra != nullptr )
                    pOutExtra->vJumpSeeds.emplace_back( *xItFirst, xDummySeed );
            } // if
            auto xItLast = pSeeds->rbegin( );
            while( xItLast != pSeeds->rend( ) && xItLast->size( ) == 0 )
                xItLast++;
            if( xItLast != pSeeds->rend( ) && xItLast->end( ) + uiMinDistDummy <= pQuery->length( ) )
            {
                pRet->emplace_back( *xItLast, pQuery->length( ), false, pQuery->iId, uiMaxDistDummy );
                if( pOutExtra != nullptr )
                    pOutExtra->vJumpSeeds.emplace_back( xDummySeed, *xItLast );
            } // if
        } // if

        return pRet;
    } // method


    std::shared_ptr<libMS::ContainerVector<SvJump>>
    computeJumpsPy( std::shared_ptr<Seeds> pSeeds, std::shared_ptr<NucSeq> pQuery, std::shared_ptr<Pack> pRefSeq )
    {
        return computeJumps( pSeeds, pQuery, pRefSeq, nullptr );
    } // method

    std::shared_ptr<Seeds> extractSeeds( std::shared_ptr<SegmentVector> pSegments, std::shared_ptr<FMIndex> pFM_index,
                                         std::shared_ptr<NucSeq> pQuery, HelperRetVal* pOutExtra )
    {
        auto pSeeds = std::make_shared<Seeds>( );
        // avoid multiple allocations (we can only guess the actual number of seeds here)
        pSeeds->reserve( pSegments->size( ) * 2 );
        pSegments->emplaceAllEachSeeds( *pFM_index, pQuery->length( ), uiMaxAmbiguitySv, uiMinSeedSizeSV, *pSeeds,
                                        [ & ]( ) { return true; } );
        if( pOutExtra != nullptr )
            xParlindromeFilter.keepParlindromes( );
        auto pFilteredSeeds = xParlindromeFilter.execute( pSeeds );
        if( pOutExtra != nullptr )
        {
            for( auto& rSeed : *pFilteredSeeds )
            {
                pOutExtra->pSeeds->push_back( rSeed );
                pOutExtra->vSocIds.push_back( pOutExtra->uiCurrSocID );
                pOutExtra->vLayerOfSeeds.push_back( 0 );
                pOutExtra->vParlindromeSeed.push_back( false );
            } // for
            for( auto& rSeed : *xParlindromeFilter.pParlindromes )
            {
                pOutExtra->pSeeds->push_back( rSeed );
                pOutExtra->vSocIds.push_back( pOutExtra->uiCurrSocID );
                pOutExtra->vLayerOfSeeds.push_back( 0 );
                pOutExtra->vParlindromeSeed.push_back( true );
            } // for
        } // if
        return pFilteredSeeds;
    }

    /**
     * @brief computes area between two seeds
     * @details
     * 1) For reseeding between two seeds we need the appropriate interval on query and reference.
     * 2) For reseeding before/after a seed we need the appropiate rectangle as well.
     *      The width of this rectangle will be its height * dExtraSeedingAreaFactor (@todo move to settings).
     *      This case is triggered by passing xDummySeed as rLast or rNext.
     *
     * If the rectangle's width is more than xMaxSizeReseed (see settings) this will return
     * two rectangles using case 2; one for each seed.
     */
    std::pair<geom::Rectangle<nucSeqIndex>, geom::Rectangle<nucSeqIndex>> DLL_PORT( MSV )
        getPositionsForSeeds( Seed& rLast, Seed& rNext, nucSeqIndex uiQStart, nucSeqIndex uiQEnd,
                              std::shared_ptr<Pack> pRefSeq );

    /**
     * @brief computes how much percent of the rectangles xRects is filled by seeds in pvSeeds.
     * @details
     * Assumes that the seeds are completeley within the rectangles.
     */
    float rectFillPercentage(
        std::shared_ptr<Seeds> pvSeeds, std::pair<geom::Rectangle<nucSeqIndex>, geom::Rectangle<nucSeqIndex>> xRects )
    {
        nucSeqIndex uiSeedSize = 0;
        for( auto& rSeed : *pvSeeds )
            uiSeedSize += rSeed.size( );
        if( xRects.first.xXAxis.size( ) * xRects.first.xYAxis.size( ) +
                xRects.second.xXAxis.size( ) * xRects.second.xYAxis.size( ) ==
            0 )
            return 0;
        return uiSeedSize / (float)( xRects.first.xXAxis.size( ) * xRects.first.xYAxis.size( ) +
                                     xRects.second.xXAxis.size( ) * xRects.second.xYAxis.size( ) );
    } // method
    float rectFillPercentage( std::shared_ptr<Seeds> pvSeeds, geom::Rectangle<nucSeqIndex> xRect )
    {
        nucSeqIndex uiSeedSize = 0;
        for( auto& rSeed : *pvSeeds )
            uiSeedSize += rSeed.size( );
        if( xRect.xXAxis.size( ) * xRect.xYAxis.size( ) == 0 )
            return 0;
        return uiSeedSize / (float)( xRect.xXAxis.size( ) * xRect.xYAxis.size( ) );
    } // method

    /**
     * @brief computes all seeds within xArea.
     * @details
     * computes all SvJumpsFromSeeds::getKMerSizeForRectangle( xArea ) sized seeds within xArea and appends
     * them to rvRet.
     * @note This is a helper function. Use the other computeSeeds.
     */
    void DLL_PORT( MSV )
        computeSeeds( geom::Rectangle<nucSeqIndex>& xArea, std::shared_ptr<NucSeq> pQuery,
                      std::shared_ptr<Pack> pRefSeq, std::shared_ptr<Seeds> rvRet, HelperRetVal* pOutExtra );
    /**
     * @brief computes all seeds within the given areas.
     * @details
     * computes all seeds larger equal to SvJumpsFromSeeds::getKMerSizeForRectangle( xArea ) within xAreas.first and
     * xAreas.second seperately.
     */
    std::shared_ptr<Seeds> DLL_PORT( MSV )
        computeSeeds( std::pair<geom::Rectangle<nucSeqIndex>, geom::Rectangle<nucSeqIndex>>& xAreas,
                      std::shared_ptr<NucSeq> pQuery, std::shared_ptr<Pack> pRefSeq, HelperRetVal* pOutExtra );

    std::shared_ptr<Seeds> DLL_PORT( MSV )
        computeSeeds( geom::Rectangle<nucSeqIndex> xArea, std::shared_ptr<NucSeq> pQuery, std::shared_ptr<Pack> pRefSeq,
                      HelperRetVal* pOutExtra );


    /**
     * @brief computes all SV jumps between the given seeds.
     * @details
     * Filters the initial seeds by their distance to the next unique seed on query:
     *      Keeps only the seeds thats closest (on the reference) to the next/last unique seed on the query.
     *
     * Uses makeJumpsByReseedingRecursive().
     */
    std::shared_ptr<libMS::ContainerVector<SvJump>> DLL_PORT( MSV )
        execute_helper( std::shared_ptr<SegmentVector> pSegments,
                        std::shared_ptr<Pack>
                            pRefSeq,
                        std::shared_ptr<FMIndex>
                            pFM_index,
                        std::shared_ptr<NucSeq>
                            pQuery,
                        HelperRetVal* pOutExtra )
    {
        return computeJumps(
            reseed( extractSeeds( pSegments, pFM_index, pQuery, pOutExtra ), pQuery, pRefSeq, pOutExtra ), pQuery,
            pRefSeq, pOutExtra );
    }

    inline HelperRetVal execute_helper_py( std::shared_ptr<SegmentVector> pSegments,
                                           std::shared_ptr<Pack>
                                               pRefSeq,
                                           std::shared_ptr<FMIndex>
                                               pFM_index,
                                           std::shared_ptr<NucSeq>
                                               pQuery )
    {
        HelperRetVal xRet;
        execute_helper( pSegments, pRefSeq, pFM_index, pQuery, &xRet );
        return xRet;
    } // method
    inline std::pair<HelperRetVal, std::shared_ptr<Seeds>>
    execute_helper_py2( std::shared_ptr<Seeds> pSeeds, std::shared_ptr<Pack> pRefSeq, std::shared_ptr<NucSeq> pQuery )
    {
        HelperRetVal xRet;
        auto pReseeded = reseed( pSeeds, pQuery, pRefSeq, &xRet );
        computeJumps( pReseeded, pQuery, pRefSeq, &xRet );
        return std::make_pair( xRet, pReseeded );
    } // method

    inline HelperRetVal execute_helper_py3( std::shared_ptr<Seeds> pSeeds, std::shared_ptr<Pack> pRefSeq,
                                            std::shared_ptr<NucSeq> pQuery )
    {
        HelperRetVal xRet;
        computeJumps( pSeeds, pQuery, pRefSeq, &xRet );
        return xRet;
    } // method

    virtual std::shared_ptr<libMS::ContainerVector<SvJump>> DLL_PORT( MSV )
        execute( std::shared_ptr<SegmentVector> pSegments,
                 std::shared_ptr<Pack>
                     pRefSeq,
                 std::shared_ptr<FMIndex>
                     pFM_index,
                 std::shared_ptr<NucSeq>
                     pQuery )
    {
        return execute_helper( pSegments, pRefSeq, pFM_index, pQuery, nullptr );
    }
}; // class


class SvJumpsFromExtractedSeeds : public libMS::Module<libMS::ContainerVector<SvJump>, false, Seeds, Pack, NucSeq>
{
    SvJumpsFromSeeds xJumpsFromSeeds;

  public:
    SvJumpsFromExtractedSeeds( const ParameterSetManager& rParameters, std::shared_ptr<Pack> pRefSeq )
        : xJumpsFromSeeds( rParameters, pRefSeq )
    {}

    virtual std::shared_ptr<libMS::ContainerVector<SvJump>>
    execute( std::shared_ptr<Seeds> pSeeds, std::shared_ptr<Pack> pRefSeq, std::shared_ptr<NucSeq> pQuery )
    {
        return xJumpsFromSeeds.computeJumps( pSeeds, pQuery, pRefSeq, nullptr );
        // return xJumpsFromSeeds.computeJumps( xJumpsFromSeeds.reseed( pSeeds, pQuery, pRefSeq, nullptr ), pQuery,
        //                                     pRefSeq, nullptr );
    }
}; // class

class ExtractSeedsFilter : public libMS::Module<Seeds, false, SegmentVector, Pack, FMIndex, NucSeq>
{
    SvJumpsFromSeeds xJumpsFromSeeds;
    nucSeqIndex uiFrom;
    nucSeqIndex uiTo;


  public:
    ExtractSeedsFilter( const ParameterSetManager& rParameters, std::shared_ptr<Pack> pRefSeq, nucSeqIndex uiFrom,
                        nucSeqIndex uiTo )
        : xJumpsFromSeeds( rParameters, pRefSeq ), uiFrom( uiFrom ), uiTo( uiTo )
    {}

    virtual std::shared_ptr<Seeds> execute( std::shared_ptr<SegmentVector> pSegments,
                                            std::shared_ptr<Pack>
                                                pRefSeq,
                                            std::shared_ptr<FMIndex>
                                                pFM_index,
                                            std::shared_ptr<NucSeq>
                                                pQuery )
    {
        auto pSeeds = xJumpsFromSeeds.extractSeeds( pSegments, pFM_index, pQuery, nullptr );

        pSeeds->erase( std::remove_if( pSeeds->begin( ), pSeeds->end( ),
                                       [ & ]( Seed& rS ) { return rS.start_ref( ) > uiTo || rS.end_ref( ) < uiFrom; } ),
                       pSeeds->end( ) );
        return pSeeds;
    }
}; // class

class RecursiveReseedingSegments : public libMS::Module<Seeds, false, SegmentVector, Pack, FMIndex, NucSeq>
{
    SvJumpsFromSeeds xJumpsFromSeeds;

  public:
    RecursiveReseedingSegments( const ParameterSetManager& rParameters, std::shared_ptr<Pack> pRefSeq )
        : xJumpsFromSeeds( rParameters, pRefSeq )
    {}

    virtual std::shared_ptr<Seeds> execute( std::shared_ptr<SegmentVector> pSegments,
                                            std::shared_ptr<Pack>
                                                pRefSeq,
                                            std::shared_ptr<FMIndex>
                                                pFM_index,
                                            std::shared_ptr<NucSeq>
                                                pQuery )
    {
        return xJumpsFromSeeds.reseed( xJumpsFromSeeds.extractSeeds( pSegments, pFM_index, pQuery, nullptr ), pQuery,
                                       pRefSeq, nullptr );
    }
}; // class

class RecursiveReseeding : public libMS::Module<Seeds, false, Seeds, Pack, NucSeq>
{
    SvJumpsFromSeeds xJumpsFromSeeds;

  public:
    RecursiveReseeding( const ParameterSetManager& rParameters, std::shared_ptr<Pack> pRefSeq )
        : xJumpsFromSeeds( rParameters, pRefSeq )
    {}

    virtual std::shared_ptr<Seeds> execute( std::shared_ptr<Seeds> pSeeds, std::shared_ptr<Pack> pRefSeq,
                                            std::shared_ptr<NucSeq> pQuery )
    {
        return xJumpsFromSeeds.reseed( pSeeds, pQuery, pRefSeq, nullptr );
    }
}; // class

class RecursiveReseedingSoCs : public libMS::Module<Seeds, false, SeedsSet, Pack, NucSeq>
{
    SvJumpsFromSeeds xJumpsFromSeeds;
    SortRemoveDuplicates xDupRem;
    FilterOverlappingSeeds xFilter;
    nucSeqIndex uiMinNt;
    const size_t uiSoCHeight;

  public:
    RecursiveReseedingSoCs( const ParameterSetManager& rParameters, std::shared_ptr<Pack> pRefSeq )
        : xJumpsFromSeeds( rParameters, pRefSeq ),
          xDupRem( rParameters ),
          xFilter( rParameters ),
          uiMinNt( rParameters.getSelected( )->xMinNtAfterReseeding->get( ) ),
          uiSoCHeight( rParameters.getSelected( )->xSoCWidth->get( ) ) // same as width
    {}

    /**
     * @note not threadsave if pOutExtra != nullptr
     */
    virtual std::shared_ptr<Seeds> execute_helper( std::shared_ptr<SeedsSet> pSeedsSet, std::shared_ptr<Pack> pRefSeq,
                                                   std::shared_ptr<NucSeq> pQuery, HelperRetVal* pOutExtra )
    {
        std::vector<std::shared_ptr<Seeds>> vSoCs;
        vSoCs.push_back( std::make_shared<Seeds>( ) );
        nucSeqIndex uiNumNtLast = 0;
        nucSeqIndex uiSizeLastSoC = 0;
        // lambda function used later (used twice avoiding duplicate code)
        auto fSoCInsertedPost = [ & ]( ) {
            // last SoC had to little NT in it
            if( uiNumNtLast < uiMinNt || uiNumNtLast == 0 )
            {
                if( pOutExtra != nullptr )
                    pOutExtra->pRemovedSeeds->append( vSoCs.back( ) );
                vSoCs.pop_back( );
            } // if
            vSoCs.push_back( std::make_shared<Seeds>( ) );
            uiNumNtLast = 0;
            uiSizeLastSoC = 0;
        };
        for( auto pSeeds : pSeedsSet->xContent )
        {
            if( pOutExtra != nullptr )
                pOutExtra->uiCurrSocID += 1;
            auto pR = xJumpsFromSeeds.reseed( pSeeds, pQuery, pRefSeq, pOutExtra );
            std::sort( pR->begin( ), pR->end( ),
                       []( const Seed& rA, const Seed& rB ) { return rA.start( ) < rB.start( ); } );
            nucSeqIndex uiMaxQ = 0;
            uiNumNtLast = 0;
            uiSizeLastSoC = 0;
            for( Seed& xSeed : *pR )
            {
                // break the SoC since there is a gap in it (seeds too far from each other on query)
                if( xSeed.start( ) > uiMaxQ + uiSoCHeight )
                    fSoCInsertedPost( );
                uiNumNtLast += xSeed.size( );
                uiMaxQ = std::max( uiMaxQ, xSeed.end( ) );
                vSoCs.back( )->push_back( xSeed );
                uiSizeLastSoC++;
            } // for
            fSoCInsertedPost( );
        } // for
        vSoCs.pop_back( );
        if( pOutExtra != nullptr )
            pOutExtra->uiCurrSocID = 0;
        std::sort( vSoCs.begin( ), vSoCs.end( ),
                   []( const auto& pA, const auto& pB ) { return pA->front( ).start( ) < pB->front( ).start( ); } );

        // remove all enclosed SoCs
        // nucSeqIndex uiMaxEnd = 0;
        auto pIntermediate = std::make_shared<Seeds>( );
        for( auto& pSoC : vSoCs )
        {
            // if( pSoC->back( ).end( ) > uiMaxEnd )
            //{
            for( auto& rSeed : *pSoC )
                pIntermediate->push_back( rSeed );
            //    uiMaxEnd = pSoC->back( ).end( );
            //} // if
            // else if( pOutExtra != nullptr )
            //    pOutExtra->pRemovedSeeds->append( pSoC );
        } // for
        // return pRet;
        auto pFiltered = xFilter.execute_helper( xDupRem.execute( pIntermediate ), pOutExtra );

        // extra step: separate by

        auto pRet = std::make_shared<Seeds>( );

        std::sort( pFiltered->begin( ), pFiltered->end( ),
                   []( auto& xA, auto& xB ) { return xA.start( ) < xB.start( ); } );
        uiNumNtLast = 0;
        nucSeqIndex uiMaxQ = 0;
        uiSizeLastSoC = 0;
        for( Seed& xSeed : *pFiltered )
        {
            if( xSeed.start( ) > uiMaxQ + uiSoCHeight )
            {
                // last SoC had to little NT in it
                if( uiNumNtLast < uiMinNt )
                    pRet->resize( pRet->size( ) - uiSizeLastSoC ); // remove it
                uiNumNtLast = 0;
                uiSizeLastSoC = 0;
            } // if
            uiNumNtLast += xSeed.size( );
            uiMaxQ = std::max( uiMaxQ, xSeed.end( ) );
            pRet->push_back( xSeed );
            uiSizeLastSoC++;
        } // for
        // very last SoC had to little NT in it
        if( uiNumNtLast < uiMinNt )
            pRet->resize( pRet->size( ) - uiSizeLastSoC ); // remove it

        if( pOutExtra != nullptr )
            pOutExtra->computeOverlappingSeedsVector( pRet );

        return pRet;
    }
    inline std::pair<std::shared_ptr<Seeds>, HelperRetVal> execute_helper_py(
        std::shared_ptr<SeedsSet> pSeedsSet, std::shared_ptr<Pack> pRefSeq, std::shared_ptr<NucSeq> pQuery )
    {
        HelperRetVal xRet;
        return std::make_pair( execute_helper( pSeedsSet, pRefSeq, pQuery, &xRet ), xRet );
    } // method

    virtual std::shared_ptr<Seeds> execute( std::shared_ptr<SeedsSet> pSeedsSet, std::shared_ptr<Pack> pRefSeq,
                                            std::shared_ptr<NucSeq> pQuery )
    {
        return execute_helper( pSeedsSet, pRefSeq, pQuery, nullptr );
    } // method
}; // class

class FilterJumpsByRegion : public libMS::Module<libMS::ContainerVector<SvJump>, false, libMS::ContainerVector<SvJump>>
{
  public:
    int64_t iFrom;
    int64_t iTo;

    FilterJumpsByRegion( const ParameterSetManager& rParameters, int64_t iFrom, int64_t iTo )
        : iFrom( iFrom ), iTo( iTo )
    {} // constructor

    std::shared_ptr<libMS::ContainerVector<SvJump>> execute( std::shared_ptr<libMS::ContainerVector<SvJump>> pJumps )
    {
        pJumps->erase( std::remove_if( pJumps->begin( ), pJumps->end( ),
                                       [ & ]( SvJump& rJ ) {
                                           if( rJ.from_start_same_strand( ) < iTo &&
                                               rJ.from_end_same_strand( ) >= iFrom )
                                               return false;
                                           if( rJ.to_start( ) < iTo && rJ.to_end( ) >= iFrom )
                                               return false;
                                           return true;
                                       } ),
                       pJumps->end( ) );
        return pJumps;
    } // method

}; // class

class FilterJumpsByRegionSquare
    : public libMS::Module<libMS::ContainerVector<SvJump>, false, libMS::ContainerVector<SvJump>>
{
  public:
    int64_t iFrom;
    int64_t iTo;

    FilterJumpsByRegionSquare( const ParameterSetManager& rParameters, int64_t iFrom, int64_t iTo )
        : iFrom( iFrom ), iTo( iTo )
    {} // constructor

    std::shared_ptr<libMS::ContainerVector<SvJump>> execute( std::shared_ptr<libMS::ContainerVector<SvJump>> pJumps )
    {
        pJumps->erase( std::remove_if( pJumps->begin( ), pJumps->end( ),
                                       [ & ]( SvJump& rJ ) {
                                           if( rJ.from_start_same_strand( ) < iTo &&
                                               rJ.from_end_same_strand( ) >= iFrom && rJ.to_start( ) < iTo &&
                                               rJ.to_end( ) >= iFrom )
                                               return false;
                                           return true;
                                       } ),
                       pJumps->end( ) );
        return pJumps;
    } // method

}; // class

class FilterJumpsByRefAmbiguity
    : public libMS::Module<libMS::ContainerVector<SvJump>, false, libMS::ContainerVector<SvJump>, Pack>,
      public AbstractFilter
{
    nucSeqIndex uiDistance;
    nucSeqIndex uiMaxRefAmbiguity;

  public:
    FilterJumpsByRefAmbiguity( const ParameterSetManager& rParameters )
        : AbstractFilter( "FilterJumpsByRefAmbiguity" ),
          uiDistance( rParameters.getSelected( )->xMaxCallSizeShortCallFilter->get( ) ),
          uiMaxRefAmbiguity( rParameters.getSelected( )->xMaxRefAmbiguityJump->get( ) )
    {} // constructor

    std::shared_ptr<libMS::ContainerVector<SvJump>>
    execute( std::shared_ptr<libMS::ContainerVector<SvJump>> pJumps, std::shared_ptr<Pack> pPack )
    {
#if ANALYZE_FILTERS
        auto uiSizeBefore = pJumps->size( );
#endif
        pJumps->erase( std::remove_if( pJumps->begin( ), pJumps->end( ),
                                       [ & ]( SvJump& rJ ) {
                                           return rJ.referenceAmbiguity( uiDistance, 5, pPack ) > uiMaxRefAmbiguity;
                                       } ),
                       pJumps->end( ) );
#if ANALYZE_FILTERS
        std::lock_guard<std::mutex> xGuard( xLock );
        uiFilterTotal += uiSizeBefore;
        uiFilterKept += pJumps->size( );
#endif
        return pJumps;
    } // method
}; // class

}; // namespace libMSV

#ifdef WITH_PYTHON
/**
 * @brief exports the SvJumpsFromSeeds @ref libMSV::Module "module" to python.
 * @ingroup export
 */
void exportSvJumpsFromSeeds( libMS::SubmoduleOrganizer& xOrganizer );
#endif