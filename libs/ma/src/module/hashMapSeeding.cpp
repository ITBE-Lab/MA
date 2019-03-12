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


std::shared_ptr<Seeds> SeedLumping::execute( std::shared_ptr<Seeds> pSeeds )
{
    // forward-strand, delta, pos on ref, pos on query, open/close seed
    std::vector<std::tuple<bool, int64_t, nucSeqIndex, bool>> vPoints;
    // order of content of tuple is importent for the sort
    static_assert( true > false );
    // true = 1 > false = 0 -> therefore open seed = true so that it is ordered first

    // reserve required memory at once
    vPoints.reserve( pSeeds->size( ) * 2 );

    // fill vPoints
    for( Seed& rSeed : *pSeeds )
    {
        vPoints.emplace_back( rSeed.bOnForwStrand, rSeed.uiPosOnReference - (int64_t)rSeed.start( ), rSeed.start_ref( ),
                              true );
        vPoints.emplace_back( rSeed.bOnForwStrand, rSeed.uiPosOnReference - (int64_t)rSeed.start( ), rSeed.end_ref( ),
                              false );
    } // for

    // sort the seeds points according to: strand, delta position, ref pos, open/close.
    std::sort( vPoints.begin( ), vPoints.end( ) );

    nucSeqIndex uiStartR = 0;
    DEBUG( int64_t iLastDelta = 0; )
    nucSeqIndex uiNumOpen = 0;
    auto pRet = std::make_shared<Seeds>( );
    for( auto& tSeed : vPoints )
    {
        if( std::get<3>( tSeed ) )
        {
            uiNumOpen++;
            if( uiNumOpen == 1 )
            {
                DEBUG( iLastDelta = std::get<1>( tSeed ); )
                uiStartR = std::get<2>( tSeed );
            } // if
            DEBUG( assert( iLastDelta == std::get<1>( tSeed ) ); )
        } // if
        else
        {
            assert( uiNumOpen > 0 );
            uiNumOpen--;
            if( uiNumOpen = 0 )
                pRet->emplace_back( uiStartR - std::get<1>( tSeed ), std::get<2>( tSeed ) - uiStartR, uiStartR,
                                    std::get<3>( tSeed ) );
            DEBUG( assert( iLastDelta == std::get<1>( tSeed ) ); )
        } // else
    } // for
    assert( uiNumOpen == 0 );

    return pRet;
} // function

#ifdef WITH_PYTHON
void exportHashMapSeeding( py::module& rxPyModuleId )
{
    // export the HashMapSeeding class
    exportModule<HashMapSeeding>( rxPyModuleId, "HashMapSeeding" );
    // export the SeedLumping class
    exportModule<SeedLumping>( rxPyModuleId, "SeedLumping" );
} // function
#endif
