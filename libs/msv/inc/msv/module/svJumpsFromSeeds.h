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


    /// @brief this class exists mereley to expose the return value of execute_helper_py to python
    class HelperRetVal
    {
      public:
        std::shared_ptr<Seeds> pSeeds;
        std::vector<size_t> vLayerOfSeeds;
        std::vector<bool> vParlindromeSeed;
        std::vector<geom::Rectangle<nucSeqIndex>> vRectangles;
        std::vector<double> vRectangleFillPercentage;
        std::vector<size_t> vRectangleReferenceAmbiguity;
        std::vector<bool> vRectangleUsedDp;
        std::vector<std::pair<Seed, Seed>> vJumpSeeds;

        HelperRetVal( ) : pSeeds( std::make_shared<Seeds>( ) ){};
    }; // class

    void forMatchingSeeds( std::shared_ptr<Seeds> pSeeds, std::function<void( Seed&, Seed& )> fOut )
    {
        for( size_t uiI = 0; uiI < pSeeds->size( ); uiI++ )
        {
            if( ( *pSeeds )[ uiI ].size( ) == 0 ) // currently at duplicate seed.
                continue;

            size_t uiJ = uiI + 1;

            // skip seeds that overlap with uiI more the 95% on query
            // (size 0 seeds are skipped automatically)
            while( uiJ < pSeeds->size( ) && overlap( ( *pSeeds )[ uiI ], ( *pSeeds )[ uiJ ] ) > 0.95 )
                uiJ++;

            // uiJ now points to the first seeds that follows uiI and is not overlapping
            size_t uiK = uiJ;
            // create entry from uiI to all seeds overlapping more than 95% with uiJ (incl. itself)
            while( uiK < pSeeds->size( ) && overlap( ( *pSeeds )[ uiK ], ( *pSeeds )[ uiJ ] ) > 0.95 )
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
        auto xItFirst = pSeeds->begin( );
        while( xItFirst != pSeeds->end( ) && xItFirst->size( ) == 0 )
            xItFirst++;
        if( xItFirst != pSeeds->end( ) )
            fOut( getPositionsForSeeds( xDummySeed, *xItFirst, 0, pQuery->length( ), pRefSeq ) );
        auto xItLast = pSeeds->rbegin( );
        while( xItLast != pSeeds->rend( ) && xItLast->size( ) == 0 )
            xItLast++;
        if( xItLast != pSeeds->rend( ) )
            fOut( getPositionsForSeeds( *xItLast, xDummySeed, 0, pQuery->length( ), pRefSeq ) );

        forMatchingSeeds( pSeeds, [&]( Seed& rA, Seed& rB ) {
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
            std::remove_if( pSeeds->begin( ), pSeeds->end( ), [&]( const Seed& xS ) { return xS.size( ) == 0; } ),
            pSeeds->end( ) );
    } // method


struct SeedCmp {
    bool operator()( const Seed& rA, const Seed& rB ) const {
                if( rA.start( ) != rB.start( ) )
                    return rA.start( ) < rB.start( );
                if( rA.start_ref( ) != rB.start_ref( ) )
                    return rA.start_ref( ) < rB.start_ref( );
                if( rA.bOnForwStrand != rB.bOnForwStrand )
                    return rA.bOnForwStrand;
                return rA.size( ) < rB.size( );
            }
};
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
        std::shared_ptr<std::set<Seed, SeedCmp>> pOutExtraDup;
        if( pOutExtra != nullptr )
        {
            pOutExtraDup = std::make_shared<std::set<Seed, SeedCmp>>( );
            for( auto& rSeed : *pSeeds )
            {
                if( pOutExtraDup->count( rSeed ) > 0 )
                    continue;
                pOutExtraDup->insert( rSeed );
                pOutExtra->pSeeds->push_back( rSeed );
                pOutExtra->vLayerOfSeeds.push_back( 0 );
                pOutExtra->vParlindromeSeed.push_back( false );
            } // for
        } // if

        size_t uiLayer = 1;
        while( uiLayer <= 100 ) // @todo emergency limit -> adjust seed sizes and rectangle size...
        {
            // sort seeds and remove duplicates (from reseeding)
            eraseMarked( pRet );
            std::sort( pRet->begin( ), pRet->end( ), []( const Seed& rA, const Seed& rB ) {
                if( rA.start( ) != rB.start( ) )
                    return rA.start( ) < rB.start( );
                if( rA.start_ref( ) != rB.start_ref( ) )
                    return rA.start_ref( ) < rB.start_ref( );
                if( rA.bOnForwStrand != rB.bOnForwStrand )
                    return rA.bOnForwStrand;
                return rA.size( ) < rB.size( );
            } );
            markDuplicates( pRet );
#if SV_JUMP_FROM_SEED_DEBUG_PRINT
            for( auto& xSeed : *pRet )
                std::cout << "reseed seed " << xSeed.start( ) << ", " << xSeed.start_ref( ) << ", " << xSeed.size( )
                          << ( xSeed.bOnForwStrand ? " forw" : " rev" ) << std::endl;
#endif
            // compute new reseeding rectangles
            std::vector<geom::Rectangle<nucSeqIndex>> vNewRects;
            collectRectangles( pRet, pQuery, pRefSeq, pOutExtra, [&]( RectPair xNew ) {
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
                        if( pOutExtraDup->count( rSeed ) > 0 )
                            continue;
                        pOutExtraDup->insert( rSeed );
                        pOutExtra->pSeeds->push_back( rSeed );
                        pOutExtra->vLayerOfSeeds.push_back( uiLayer );
                        pOutExtra->vParlindromeSeed.push_back( false );
                    } // for
                    for( auto& rSeed : *xParlindromeFilter.pParlindromes )
                    {
                        if( pOutExtraDup->count( rSeed ) > 0 )
                            continue;
                        pOutExtraDup->insert( rSeed );
                        pOutExtra->pSeeds->push_back( rSeed );
                        pOutExtra->vLayerOfSeeds.push_back( uiLayer );
                        pOutExtra->vParlindromeSeed.push_back( true );
                    } // for
                } // if
            } // for
#if SV_JUMP_FROM_SEED_DEBUG_PRINT
            std::cout << " --- Layer --- " << std::endl;
#endif
            uiLayer++;
        } // while

        // std::cout << uiLayer << " " << pRet->size( ) << std::endl;
        return xFilterOverlapping.execute( xToSMEM.execute( pRet ) );
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
        forMatchingSeeds( pSeeds, [&]( Seed& rA, Seed& rB ) {
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
                                        [&]( ) { return true; } );
        if( pOutExtra != nullptr )
            xParlindromeFilter.keepParlindromes( );
        auto pFilteredSeeds = xParlindromeFilter.execute( pSeeds );
        if( pOutExtra != nullptr )
        {
            for( auto& rSeed : *pFilteredSeeds )
            {
                pOutExtra->pSeeds->push_back( rSeed );
                pOutExtra->vLayerOfSeeds.push_back( 0 );
                pOutExtra->vParlindromeSeed.push_back( false );
            } // for
            for( auto& rSeed : *xParlindromeFilter.pParlindromes )
            {
                pOutExtra->pSeeds->push_back( rSeed );
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
    float rectFillPercentage( std::shared_ptr<Seeds> pvSeeds,
                              std::pair<geom::Rectangle<nucSeqIndex>, geom::Rectangle<nucSeqIndex>> xRects )
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
    inline HelperRetVal execute_helper_py2( std::shared_ptr<Seeds> pSeeds, std::shared_ptr<Pack> pRefSeq,
                                            std::shared_ptr<NucSeq> pQuery )
    {
        HelperRetVal xRet;
        computeJumps( reseed( pSeeds, pQuery, pRefSeq, &xRet ), pQuery, pRefSeq, &xRet );
        return xRet;
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
        return xJumpsFromSeeds.computeJumps( xJumpsFromSeeds.reseed( pSeeds, pQuery, pRefSeq, nullptr ), pQuery,
                                             pRefSeq, nullptr );
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
                                       [&]( Seed& rS ) { return rS.start_ref( ) > uiTo || rS.end_ref( ) < uiFrom; } ),
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
                                       [&]( SvJump& rJ ) {
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
                                       [&]( SvJump& rJ ) {
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

    std::shared_ptr<libMS::ContainerVector<SvJump>> execute( std::shared_ptr<libMS::ContainerVector<SvJump>> pJumps,
                                                             std::shared_ptr<Pack> pPack )
    {
#if ANALYZE_FILTERS
        auto uiSizeBefore = pJumps->size( );
#endif
        pJumps->erase( std::remove_if( pJumps->begin( ), pJumps->end( ),
                                       [&]( SvJump& rJ ) {
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