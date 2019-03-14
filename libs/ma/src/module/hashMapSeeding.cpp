/**
 * @file hash_map_seeding.cpp
 * @author Markus Schmidt
 */
#include "module/hashMapSeeding.h"
#include <limits>
#include <unordered_map>

using namespace libMA;

using namespace libMA::defaults;


std::shared_ptr<Seeds> HashMapSeeding::execute( std::shared_ptr<NucSeq> pQ1, std::shared_ptr<NucSeq> pQ2 )
{
    std::unordered_multimap<std::string, size_t> xIndex;

    // insert seeds into index
    for( size_t uiI = 0; uiI < pQ2->length( ); uiI += 1 )
        xIndex.emplace( pQ2->fromTo( uiI, uiI + this->uiSeedSize ), uiI );

    auto pSeeds = std::make_shared<Seeds>( );
    for( size_t uiI = 0; uiI < pQ1->length( ); uiI += 1 )
    {
        auto tuiRange = xIndex.equal_range( pQ1->fromTo( uiI, uiI + this->uiSeedSize ) );
        for( auto xIt = tuiRange.first; xIt != tuiRange.second; ++xIt )
            pSeeds->emplace_back( uiI, this->uiSeedSize, xIt->second, true );
    } // for

    return pSeeds;
} // function

std::shared_ptr<Seeds> ReSeeding::execute( std::shared_ptr<Seeds> pSeeds, std::shared_ptr<NucSeq> pQuery,
                                           std::shared_ptr<Pack> pPack )
{
    std::sort( pSeeds->begin( ), pSeeds->end( ),
               []( const Seed& rA, const Seed& rB ) { return rA.start_ref( ) < rB.start_ref( ); } ); // sort call

    auto pCollection = std::make_shared<Seeds>( );

    for( size_t uiI = 0; uiI < pSeeds->size( ) - 1; uiI++ )
    {
        Seed& rA = ( *pSeeds )[ uiI ];
        Seed& rB = ( *pSeeds )[ uiI + 1 ];
        assert( rA.bOnForwStrand == rB.bOnForwStrand );
        if( rA.end( ) + xHashMapSeeder.uiSeedSize <= rB.start( ) &&
            rA.end_ref( ) + xHashMapSeeder.uiSeedSize <= rB.start_ref( ) )
        {
            // re-seed in the gap between seed uiI and seed uiI + 1
            auto pAppend = xHashMapSeeder.execute(
                std::make_shared<NucSeq>(
                    rA.bOnForwStrand
                        ? pQuery->fromTo( rA.end( ), rB.start( ) )
                        : pQuery->fromToComplement( pQuery->length( ) - rB.start( ), pQuery->length( ) - rA.end( ) ) ),
                pPack->vExtract( rA.end_ref( ), rB.start_ref( ) ) );
            for( Seed& rSeed : *pAppend )
            {
                rSeed.bOnForwStrand = rA.bOnForwStrand;
                rSeed.uiPosOnReference += rA.end_ref( );
                rSeed.iStart += rA.end( );
            } // for
            pCollection->append( pAppend );
        } // if
    } // for

    // re-seed before the first seed
    Seed& rFirst = pSeeds->front( );
    if( xHashMapSeeder.uiSeedSize <= rFirst.start( ) )
    {
        auto pQ1 = std::make_shared<NucSeq>(
            rFirst.bOnForwStrand
                ? pQuery->fromTo( uiPadding < rFirst.start( ) ? rFirst.start( ) - uiPadding : 0, rFirst.start( ) )
                : pQuery->fromToComplement(
                      pQuery->length( ) - rFirst.start( ),
                      std::min( pQuery->length( ), pQuery->length( ) - ( rFirst.start( ) - uiPadding ) ) ) );
        auto pQ2 = pPack->vExtract( uiPadding < rFirst.start_ref( ) ? rFirst.start_ref( ) - uiPadding : 0,
                                    rFirst.start_ref( ) );
        assert(pQ1->length() <= uiPadding);
        assert(pQ2->length() <= uiPadding);
        auto pAppend = xHashMapSeeder.execute( pQ1, pQ2 );
        for( Seed& rSeed : *pAppend )
        {
            rSeed.bOnForwStrand = rFirst.bOnForwStrand;
            if( rFirst.bOnForwStrand && uiPadding < rFirst.start( ) )
                rSeed.iStart += rFirst.start( ) - uiPadding;
            rSeed.uiPosOnReference += uiPadding < rFirst.start_ref( ) ? rFirst.start_ref( ) - uiPadding : 0;
        } // for
        pCollection->append( pAppend );
    } // if

    // re-seed after the last seed
    Seed& rLast = pSeeds->back( );
    if( rLast.end_ref( ) + xHashMapSeeder.uiSeedSize <= pQuery->length( ) )
    {
        auto pQ1 = std::make_shared<NucSeq>(
            rLast.bOnForwStrand
                ? pQuery->fromTo( rLast.end( ), std::min( pQuery->length( ), rLast.end( ) + uiPadding ) )
                : pQuery->fromToComplement( uiPadding < pQuery->length( ) - rLast.end( )
                                                ? pQuery->length( ) - ( rLast.end( ) + uiPadding )
                                                : 0,
                                            pQuery->length( ) - rLast.end( ) ) );
        auto pQ2 = pPack->vExtract( rLast.end_ref( ),
                                    std::min( rLast.end_ref( ) + uiPadding, pPack->uiStartOfReverseStrand( ) ) );
        assert(pQ1->length() <= uiPadding);
        assert(pQ2->length() <= uiPadding);
        auto pAppend = xHashMapSeeder.execute( pQ1, pQ2 );
        for( Seed& rSeed : *pAppend )
        {
            rSeed.bOnForwStrand = rLast.bOnForwStrand;
            rSeed.uiPosOnReference += rLast.end_ref( );
            rSeed.iStart += rLast.end( );
        } // for
        pCollection->append( pAppend );
    } // if

    // append the original seeds
    pCollection->append( pSeeds );

    // lump everything
    return xLumper.execute( pCollection );
} // function

#ifdef WITH_PYTHON
void exportHashMapSeeding( py::module& rxPyModuleId )
{
    // export the HashMapSeeding class
    exportModule<HashMapSeeding>( rxPyModuleId, "HashMapSeeding" );
    // export the ReSeeding class
    exportModule<ReSeeding>( rxPyModuleId, "ReSeeding" );
    // export the FillSeedSet class
    exportModule<FillSeedSet>( rxPyModuleId, "FillSeedSet" );
    // export the ExtractFilledSeedSets class
    exportModule<ExtractFilledSeedSets>( rxPyModuleId, "ExtractFilledSeedSets" );
} // function
#endif
