
/**
 * @file renderer_speedup.cpp
 * @author Markus Schmidt
 */

#ifdef WITH_PYTHON
#include "ma/container/minimizer_index.h"
#include "ma/container/pack.h"
#include "ma/module/seedFilters.h"
#include "ma/module/stripOfConsideration.h"
#include "ms/container/sv_db/pool_container.h"
#include "ms/util/pybind11.h"
#include "msv/container/sv_db/tables/kMerFilter.h"
#include "msv/container/sv_db/tables/read.h"
#include "msv/module/count_k_mers.h"
#include "msv/module/svJumpsFromSeeds.h"

using namespace libMSV;
using namespace libMS;

struct SeedInfo
{
    float fCenter;
    int64_t iReadId;
    std::string sReadName;
    nucSeqIndex uiSize;
    nucSeqIndex uiQ;
    nucSeqIndex uiR;
    nucSeqIndex uiL;
    size_t uiSeedOrderOnQuery;
    bool bOnForward;
    int64_t uiLayer;
    bool bParlindrome;
    std::pair<nucSeqIndex, nucSeqIndex> xX;
    std::pair<nucSeqIndex, nucSeqIndex> xY;
    size_t uiCategory;
}; // struct

struct RectangleInfo
{
    std::vector<geom::Rectangle<nucSeqIndex>> vRectangles;
    std::vector<double> vRectangleFillPercentage;
    std::vector<size_t> vRectangleReferenceAmbiguity;
    std::vector<bool> vRectangleUsedDp;
    size_t uiCategory;
    size_t uiEndColumnSize;
    int64_t iReadId;
}; // struct


void addSeed( Seed& rSeed,
              std::vector<SeedInfo>& vRet,
              std::vector<std::pair<nucSeqIndex, int64_t>>& vEndColumn,
              std::vector<int64_t>
                  vAllColIds,
              size_t uiCategoryCounter,
              bool bParlindrome,
              int64_t iLayer,
              int64_t iReadId,
              size_t uiSeedOrderOnQuery,
              std::string sReadName )
{
    float fCenter = rSeed.start_ref( ) + rSeed.size( ) / 2.0f;
    nucSeqIndex uiR = rSeed.start_ref( );
    auto xX = std::make_pair( rSeed.start_ref( ), rSeed.start_ref( ) + rSeed.size( ) );
    if( !rSeed.bOnForwStrand )
    {
        uiR = rSeed.start_ref( ) - rSeed.size( ) + 1;
        fCenter = rSeed.start_ref( ) - rSeed.size( ) / 2.0f + 1;
        xX = std::make_pair( rSeed.start_ref( ) + 1, rSeed.start_ref( ) - rSeed.size( ) + 1 );
    } // if

    // picks y-axis column to render seed in
    size_t uiCurrColumn = 0;
    while( uiCurrColumn < vEndColumn.size( ) &&
           uiR <= std::get<0>( vEndColumn[ uiCurrColumn ] ) +
                      ( std::get<1>( vEndColumn[ uiCurrColumn ] ) == iReadId ? 0 : 100 ) )
        uiCurrColumn++;
    if( uiCurrColumn >= vEndColumn.size( ) )
    {
        vEndColumn.emplace_back( 0, 0 );
        vAllColIds.push_back( uiCurrColumn + uiCategoryCounter );
    } // if
    std::get<0>( vEndColumn[ uiCurrColumn ] ) = uiR;
    std::get<1>( vEndColumn[ uiCurrColumn ] ) = iReadId;

    vRet.emplace_back( SeedInfo{fCenter, iReadId, sReadName, rSeed.size( ), rSeed.start( ), uiR, rSeed.size( ),
                                uiSeedOrderOnQuery, rSeed.bOnForwStrand, iLayer, bParlindrome, xX,
                                std::make_pair( rSeed.start( ), rSeed.end( ) ), uiCurrColumn + uiCategoryCounter} );
} // function

void addRectangle( std::vector<RectangleInfo>& vRectangles,
                   SvJumpsFromSeeds::HelperRetVal& xHelper,
                   size_t uiCategoryCounter,
                   int64_t iReadId,
                   size_t uiEndColumnSize )
{
    RectangleInfo xInfo{};
    xInfo.vRectangles.swap( xHelper.vRectangles );
    xInfo.vRectangleFillPercentage.swap( xHelper.vRectangleFillPercentage );
    xInfo.vRectangleReferenceAmbiguity.swap( xHelper.vRectangleReferenceAmbiguity );
    xInfo.vRectangleUsedDp.swap( xHelper.vRectangleUsedDp );
    xInfo.uiCategory = uiCategoryCounter;
    xInfo.uiEndColumnSize = uiEndColumnSize;
    xInfo.iReadId = iReadId;
    vRectangles.push_back( xInfo );
} // function


struct ReadInfo
{
    std::vector<SeedInfo> vRet;
    std::vector<RectangleInfo> vRectangles;
    std::vector<std::pair<size_t, int64_t>> vReadsNCols;
    std::vector<std::shared_ptr<NucSeq>> vReads;
}; // struct
template <typename DBCon>
std::shared_ptr<ReadInfo> seedDisplaysForReadIds( const ParameterSetManager& rParameters,
                                                  std::shared_ptr<libMS::PoolContainer<DBCon>>
                                                      pConPool,
                                                  std::set<int64_t>
                                                      xReadIds,
                                                  std::shared_ptr<Pack>
                                                      pPack,
                                                  std::shared_ptr<minimizer::Index>
                                                      pMMIndex,
                                                  std::shared_ptr<HashCounter>
                                                      pHashCounter,
                                                  bool bDoCompress,
                                                  size_t iMaxTime )
{
    std::chrono::system_clock::time_point xEndTime =
        std::chrono::system_clock::now( ) + std::chrono::seconds( iMaxTime );

    std::vector<int64_t> vReadIds;
    vReadIds.insert( vReadIds.begin( ), xReadIds.begin( ), xReadIds.end( ) );
    std::sort( vReadIds.begin( ), vReadIds.end( ), []( int64_t iA, int64_t iB ) { return iB < iA; } );

    std::mutex xLock;
    std::vector<SeedInfo> vRet;
    std::vector<RectangleInfo> vRectangles;
    std::vector<std::pair<size_t, int64_t>> vReadsNCols;
    std::vector<std::shared_ptr<NucSeq>> vReads;
    vRet.reserve( vReadIds.size( ) * 500 );

    MMFilteredSeeding xSeeding( rParameters, 300 );
    SeedLumping xLumping( rParameters );
    StripOfConsiderationSeeds xSoc( rParameters );
    GetAllFeasibleSoCs xSocFilter( rParameters, 100 );

    std::vector<std::future<void>> vFutures;
    // seed_order_on_query, seed, layer, parlindrome, read_id, read_name
    std::vector<std::tuple<size_t, Seed, size_t, bool, int64_t, std::string>> vAllSeeds;
    // x-end position of last seeds and the seeds read id for each column
    std::vector<int64_t> vAllColIds;
    size_t uiCategoryCounter = 0;
    bool bStop = false;
    for( size_t uiI = 0; uiI < rParameters.getNumThreads( ); uiI++ )
        vFutures.push_back( pConPool->xPool.enqueue(
            [&]( std::shared_ptr<DBCon> pConn, size_t uiI ) { //
                SvJumpsFromSeeds xJumpsFromSeeds( rParameters, pPack );
                auto pReadTable = std::make_shared<ReadTable<DBCon>>( pConn );
                for( size_t uiJ = uiI; uiJ < vReadIds.size( ); uiJ += rParameters.getNumThreads( ) )
                {
                    auto iReadId = vReadIds[ uiJ ];
                    auto pRead = pReadTable->getRead( iReadId );
                    auto pMinimizers = xSeeding.execute( pMMIndex, pRead, pPack, pHashCounter );
                    auto pLumpedSeeds = xLumping.execute( pMinimizers, pRead, pPack );
                    auto pSoCs = xSoc.execute( pLumpedSeeds, pRead, pPack );
                    auto pFilteredSeeds = xSocFilter.execute( pSoCs );

                    auto xHelperRet = xJumpsFromSeeds.execute_helper_py2( pFilteredSeeds, pPack, pRead );

                    // seed_order_on_query, seed, layer, parlindrome, read_id, read_name
                    std::vector<std::tuple<size_t, Seed, size_t, bool, int64_t, std::string>> vSeedsNIndex;
                    for( size_t uiK = 0; uiK < xHelperRet.pSeeds->size( ); uiK++ )
                        vSeedsNIndex.emplace_back( 0,
                                                   ( *xHelperRet.pSeeds )[ uiK ],
                                                   xHelperRet.vLayerOfSeeds[ uiK ],
                                                   xHelperRet.vParlindromeSeed[ uiK ],
                                                   iReadId,
                                                   pRead->sName );
                    std::sort( vSeedsNIndex.begin( ), vSeedsNIndex.end( ), []( auto& xA, auto& xB ) {
                        return std::get<1>( xA ).start( ) < std::get<1>( xB ).start( );
                    } );
                    for( size_t uiK = 0; uiK < vSeedsNIndex.size( ); uiK++ )
                        std::get<0>( vSeedsNIndex[ uiK ] ) = uiK;


                    if( bDoCompress )
                    {
                        std::lock_guard<std::mutex> xGuard( xLock );
                        if( bStop )
                            return;
                        vAllSeeds.insert( vAllSeeds.end( ), vSeedsNIndex.begin( ), vSeedsNIndex.end( ) );
                        addRectangle( vRectangles, xHelperRet, uiCategoryCounter, iReadId, 0 );
                        vReads.push_back( pRead );
                    } // if
                    else
                    {
                        std::sort( vSeedsNIndex.begin( ), vSeedsNIndex.end( ), []( auto& xA, auto& xB ) {
                            return std::get<1>( xA ).start_ref( ) < std::get<1>( xB ).start_ref( );
                        } );
                        std::vector<std::pair<nucSeqIndex, int64_t>> vEndColumn;
                        std::lock_guard<std::mutex> xGuard( xLock );
                        if( bStop )
                            return;
                        for( auto& xTup : vAllSeeds )
                            addSeed( std::get<1>( xTup ),
                                     vRet,
                                     vEndColumn,
                                     vAllColIds,
                                     uiCategoryCounter,
                                     std::get<3>( xTup ),
                                     std::get<2>( xTup ),
                                     std::get<4>( xTup ),
                                     std::get<0>( xTup ),
                                     std::get<5>( xTup ) );
                        vAllColIds.push_back( uiCategoryCounter + ( vEndColumn.size( ) - 1 ) / 2 );

                        addRectangle( vRectangles, xHelperRet, uiCategoryCounter, iReadId, vEndColumn.size( ) );

                        uiCategoryCounter += vEndColumn.size( );
                        vReadsNCols.emplace_back( vAllColIds.back( ), iReadId );
                        vReads.push_back( pRead );
                    } // else
                } // for
            },
            uiI ) );

    // wait for threads to finish at most iMaxTime seconds, then stop all work and give up
    for( auto& xFuture : vFutures )
    {
        // get status until it is either timeout or ready
        std::future_status xStatus = std::future_status::deferred;
        while( xStatus == std::future_status::deferred )
            xStatus = xFuture.wait_until( xEndTime );

        if( xStatus == std::future_status::timeout )
        {
            std::lock_guard<std::mutex> xGuard( xLock );
            bStop = true;
        } // if
        else if( xStatus == std::future_status::ready )
            xFuture.get( );
        else
            throw std::runtime_error( "should be unreachable" );
    } // for

    if( bStop )
    {
        vRet.clear( );
        vRectangles.clear( );
        vReadsNCols.clear( );
        vReads.clear( );
    } // if
    else if( bDoCompress && !vAllSeeds.empty( ) )
    {
        std::vector<std::pair<nucSeqIndex, int64_t>> vEndColumn;
        std::sort( vAllSeeds.begin( ), vAllSeeds.end( ), []( auto& xA, auto& xB ) {
            if( std::get<1>( xA ).start_ref( ) != std::get<1>( xB ).start_ref( ) )
                return std::get<1>( xA ).start_ref( ) < std::get<1>( xB ).start_ref( );
            return std::get<4>( xA ) < std::get<4>( xB );
        } );
        for( auto& xTup : vAllSeeds )
            addSeed( std::get<1>( xTup ),
                     vRet,
                     vEndColumn,
                     vAllColIds,
                     uiCategoryCounter,
                     std::get<3>( xTup ),
                     std::get<2>( xTup ),
                     std::get<4>( xTup ),
                     std::get<0>( xTup ),
                     std::get<5>( xTup ) );
        uiCategoryCounter += vEndColumn.size( );
    } // if

    return std::make_shared<ReadInfo>( ReadInfo{vRet, vRectangles, vReadsNCols, vReads} );
} // method


#include "ms/container/sv_db/py_db_conf.h"
#include <pybind11/stl.h>
void exportRendererSpeedUp( libMS::SubmoduleOrganizer& xOrganizer )
{
    py::class_<SeedInfo>( xOrganizer._util( ), "SeedInfo" )
        .def_readwrite( "fCenter", &SeedInfo::fCenter )
        .def_readwrite( "iReadId", &SeedInfo::iReadId )
        .def_readwrite( "sReadName", &SeedInfo::sReadName )
        .def_readwrite( "uiSize", &SeedInfo::uiSize )
        .def_readwrite( "uiQ", &SeedInfo::uiQ )
        .def_readwrite( "uiR", &SeedInfo::uiR )
        .def_readwrite( "uiL", &SeedInfo::uiL )
        .def_readwrite( "uiSeedOrderOnQuery", &SeedInfo::uiSeedOrderOnQuery )
        .def_readwrite( "bOnForward", &SeedInfo::bOnForward )
        .def_readwrite( "uiLayer", &SeedInfo::uiLayer )
        .def_readwrite( "bParlindrome", &SeedInfo::bParlindrome )
        .def_readwrite( "xX", &SeedInfo::xX )
        .def_readwrite( "xY", &SeedInfo::xY )
        .def_readwrite( "uiCategory", &SeedInfo::uiCategory );
    py::class_<RectangleInfo>( xOrganizer._util( ), "RectangleInfo" )
        .def_readwrite( "vRectangles", &RectangleInfo::vRectangles )
        .def_readwrite( "vRectangleFillPercentage", &RectangleInfo::vRectangleFillPercentage )
        .def_readwrite( "vRectangleReferenceAmbiguity", &RectangleInfo::vRectangleReferenceAmbiguity )
        .def_readwrite( "vRectangleUsedDp", &RectangleInfo::vRectangleUsedDp )
        .def_readwrite( "uiCategory", &RectangleInfo::uiCategory )
        .def_readwrite( "iReadId", &RectangleInfo::iReadId )
        .def_readwrite( "uiEndColumnSize", &RectangleInfo::uiEndColumnSize );
    py::class_<ReadInfo, std::shared_ptr<ReadInfo>>( xOrganizer._util( ), "ReadInfo" )
        .def_readwrite( "vRet", &ReadInfo::vRet )
        .def_readwrite( "vRectangles", &ReadInfo::vRectangles )
        .def_readwrite( "vReadsNCols", &ReadInfo::vReadsNCols )
        .def_readwrite( "vReads", &ReadInfo::vReads );

    xOrganizer.util( ).def( "seedDisplaysForReadIds", &seedDisplaysForReadIds<DBCon> );
} // function
#endif