
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
    bool bOverlapping;
    bool bInSocReseed;
    std::pair<nucSeqIndex, nucSeqIndex> xX;
    std::pair<nucSeqIndex, nucSeqIndex> xY;
    size_t uiCategory;
    size_t uiMinFilterCount;
    size_t uiMaxFilterCount;
    size_t uiSoCNt;
    size_t uiSocId;
}; // struct

struct RectangleInfo
{
    std::vector<geom::Rectangle<nucSeqIndex>> vRectangles;
    std::vector<size_t> vRectangleLayers;
    std::vector<double> vRectangleFillPercentage;
    std::vector<size_t> vRectangleReferenceAmbiguity;
    std::vector<size_t> vRectangleKMerSize;
    std::vector<bool> vRectangleUsedDp;
    size_t uiCategory;
    size_t uiEndColumnSize;
    int64_t iReadId;
    bool bInSoCReseeding;
}; // struct


void addSeed( Seed& rSeed, std::vector<SeedInfo>& vRet, std::vector<std::pair<nucSeqIndex, int64_t>>& vEndColumn,
              std::vector<int64_t>& vAllColIds, size_t uiCategoryCounter, bool bParlindrome, bool bOverlapping,
              int64_t iLayer, int64_t iReadId, size_t uiSeedOrderOnQuery, std::string sReadName,
              std::shared_ptr<HashCounter> pCounter, std::shared_ptr<Pack> pPack,
              std::shared_ptr<minimizer::Index> pMMIndex, std::shared_ptr<NucSeq> pRead, bool bInSocReseed,
              bool bRectSoc, size_t uiSocId )
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
    std::get<0>( vEndColumn[ uiCurrColumn ] ) = uiR + rSeed.size( );
    std::get<1>( vEndColumn[ uiCurrColumn ] ) = iReadId;


    size_t uiMax = MMFilteredSeeding::getMinCount( rSeed, pMMIndex, pRead, pPack, pCounter, bRectSoc );
    size_t uiMin = MMFilteredSeeding::getMaxCount( rSeed, pMMIndex, pRead, pPack, pCounter, bRectSoc );

    vRet.emplace_back( SeedInfo{ fCenter,
                                 iReadId,
                                 sReadName,
                                 rSeed.size( ),
                                 rSeed.start( ),
                                 uiR,
                                 rSeed.size( ),
                                 uiSeedOrderOnQuery,
                                 rSeed.bOnForwStrand,
                                 iLayer,
                                 bParlindrome,
                                 bOverlapping,
                                 bInSocReseed,
                                 xX,
                                 std::make_pair( rSeed.start( ), rSeed.end( ) ),
                                 uiCurrColumn + uiCategoryCounter,
                                 uiMin,
                                 uiMax,
                                 rSeed.uiSoCNt,
                                 uiSocId } );
} // function

void addRectangle( std::vector<RectangleInfo>& vRectangles,
                   HelperRetVal& xHelper,
                   size_t uiCategoryCounter,
                   int64_t iReadId,
                   size_t uiEndColumnSize,
                   bool bInSoCReseeding )
{
    RectangleInfo xInfo{ };
    xInfo.vRectangles.swap( xHelper.vRectangles );
    xInfo.vRectangleLayers.swap( xHelper.vRectangleLayers );
    xInfo.vRectangleFillPercentage.swap( xHelper.vRectangleFillPercentage );
    xInfo.vRectangleReferenceAmbiguity.swap( xHelper.vRectangleReferenceAmbiguity );
    xInfo.vRectangleKMerSize.swap( xHelper.vRectangleKMerSize );
    xInfo.vRectangleUsedDp.swap( xHelper.vRectangleUsedDp );
    xInfo.uiCategory = uiCategoryCounter;
    xInfo.uiEndColumnSize = uiEndColumnSize;
    xInfo.iReadId = iReadId;
    xInfo.bInSoCReseeding = bInSoCReseeding;
    vRectangles.push_back( xInfo );
} // function


struct ReadInfo
{
    std::vector<SeedInfo> vRet;
    std::vector<RectangleInfo> vRectangles;
    std::vector<std::pair<size_t, int64_t>> vReadsNCols;
    std::vector<std::shared_ptr<NucSeq>> vReads;
    std::vector<int64_t> vColIds;
    std::vector<int64_t> vAllColIds;
}; // struct

struct HashCounters
{
    std::map<int64_t, std::shared_ptr<HashCounter>> xCounters;
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
                                                  std::shared_ptr<HashCounters>
                                                      pHashCounters,
                                                  bool bDoCompress,
                                                  size_t iMaxTime )
{
    std::chrono::system_clock::time_point xEndTime =
        std::chrono::system_clock::now( ) + std::chrono::milliseconds( iMaxTime );

    std::vector<int64_t> vReadIds;
    vReadIds.insert( vReadIds.begin( ), xReadIds.begin( ), xReadIds.end( ) );
    std::sort( vReadIds.begin( ), vReadIds.end( ), []( int64_t iA, int64_t iB ) { return iB < iA; } );

    std::mutex xLock;
    std::mutex xLock2;
    std::vector<SeedInfo> vRet;
    std::vector<RectangleInfo> vRectangles;
    std::vector<std::pair<size_t, int64_t>> vReadsNCols;
    std::vector<std::shared_ptr<NucSeq>> vReads;
    vRet.reserve( vReadIds.size( ) * 500 );

    MMFilteredSeeding xSeeding( rParameters );
    SeedLumping xLumping( rParameters );
    //FilterContigBorder xCtgFilter( rParameters );
    StripOfConsiderationSeeds xSoc( rParameters );
    GetAllFeasibleSoCsAsSet xSocFilter( rParameters );

    std::vector<std::future<void>> vFutures;
    // seed_order_on_query, seed, layer, parlindrome, overlapping, read_id, read_name, read, bInSocReseed
    std::vector<std::tuple<size_t, Seed, size_t, bool, bool, int64_t, std::string, std::shared_ptr<NucSeq>, bool,
                           std::shared_ptr<HashCounter>, size_t>>
        vAllSeeds;
    // x-end position of last seeds and the seeds read id for each column
    std::vector<int64_t> vColIds;
    std::vector<int64_t> vAllColIds;
    size_t uiCategoryCounter = 0;
    bool bStop = false;

    for( size_t uiI = 0; uiI < rParameters.getNumThreads( ); uiI++ )
        vFutures.push_back( pConPool->xPool.enqueue(
            [ & ]( std::shared_ptr<DBCon> pConn, size_t uiI ) { //
                auto pReadTable = std::make_unique<ReadTable<DBCon>>( pConn );
                for( size_t uiJ = uiI; uiJ < vReadIds.size( ); uiJ += rParameters.getNumThreads( ) )
                {
                    // modules not threadsave if HelperRetVal is given
                    RecursiveReseedingSoCs xReseeding( rParameters, pPack );
                    SvJumpsFromSeeds xJumpsFromSeeds( rParameters, pPack );
                    auto iReadId = vReadIds[ uiJ ];
                    int64_t uiSeqId = pReadTable->getSeqId( iReadId );
                    std::shared_ptr<HashCounter> pCounter;
                    {
                        std::lock_guard<std::mutex> xGuard2( xLock2 );
                        auto pIt = pHashCounters->xCounters.find( uiSeqId );
                        if( pIt != pHashCounters->xCounters.end( ) )
                            pCounter = pIt->second;
                        else
                        {
                            pCounter = std::make_unique<HashFilterTable<DBCon>>( pConn )->getCounter( uiSeqId );
                            pHashCounters->xCounters[ uiSeqId ] = pCounter;
                        } // else
                    } // scope of xGuard2
                    auto pRead = pReadTable->getRead( iReadId );
                    auto pMinimizers = xSeeding.execute( pMMIndex, pRead, pPack, pCounter );
                    auto pLumpedSeeds = xLumping.execute( pMinimizers, pRead, pPack );
                    //auto pCtgFilteredSeeds = xCtgFilter.execute( pLumpedSeeds, pPack );
                    auto pSoCs = xSoc.execute( pLumpedSeeds, pRead, pPack );
                    auto pFilteredSeeds = xSocFilter.execute( pSoCs );
                    HelperRetVal xReseedOutExtraInfo;
                    auto pReseeded = xReseeding.execute_helper( pFilteredSeeds, pPack, pRead, &xReseedOutExtraInfo );

                    auto xHelperRet = xJumpsFromSeeds.execute_helper_py3( pReseeded, pPack, pRead );

                    // seed_order_on_query, seed, layer, parlindrome, overlapping, read_id, read_name, read,
                    // bInSocReseed
                    std::vector<std::tuple<size_t, Seed, size_t, bool, bool, int64_t, std::string,
                                           std::shared_ptr<NucSeq>, bool, std::shared_ptr<HashCounter>, size_t>>
                        vSeedsNIndex;
                    //for( size_t uiK = 0; uiK < xReseedOutExtraInfo.pSeeds->size( ); uiK++ )
                    //    if( !xReseedOutExtraInfo.vOverlappingSeed[ uiK ] )
                    //        vSeedsNIndex.emplace_back(
                    //            0, ( *xReseedOutExtraInfo.pSeeds )[ uiK ], xReseedOutExtraInfo.vLayerOfSeeds[ uiK ],
                    //            xReseedOutExtraInfo.vParlindromeSeed[ uiK ],
                    //            xReseedOutExtraInfo.vOverlappingSeed[ uiK ], iReadId, pRead->sName, pRead, true,
                    //            pCounter, xReseedOutExtraInfo.vSocIds[ uiK ] );
                    for( auto& rSeed : *xReseedOutExtraInfo.pRemovedSeeds )
                        vSeedsNIndex.emplace_back( 0, rSeed, 0, false, true, iReadId, pRead->sName, pRead, true,
                                                   pCounter, 0 );
                    for( auto& rSeed : *pReseeded )
                        vSeedsNIndex.emplace_back( 0, rSeed, 0, false, false, iReadId, pRead->sName, pRead, true,
                                                   pCounter, 0 );
                    for( size_t uiK = 0; uiK < xHelperRet.pSeeds->size( ); uiK++ )
                        // only use seeds that we do not get from xReseedOutExtraInfo already
                        if( xHelperRet.vLayerOfSeeds[ uiK ] > 0 )
                            vSeedsNIndex.emplace_back(
                                0, ( *xHelperRet.pSeeds )[ uiK ], xHelperRet.vLayerOfSeeds[ uiK ],
                                xHelperRet.vParlindromeSeed[ uiK ], xHelperRet.vOverlappingSeed[ uiK ], iReadId,
                                pRead->sName, pRead, false, pCounter, xHelperRet.vSocIds[ uiK ] );

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
                        addRectangle( vRectangles, xHelperRet, uiCategoryCounter, iReadId, 0, false );
                        addRectangle( vRectangles, xReseedOutExtraInfo, uiCategoryCounter, iReadId, 0, true );
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
                        for( auto& xTup : vSeedsNIndex )
                            addSeed( std::get<1>( xTup ),
                                     vRet,
                                     vEndColumn,
                                     vAllColIds,
                                     uiCategoryCounter,
                                     std::get<3>( xTup ),
                                     std::get<4>( xTup ),
                                     std::get<2>( xTup ),
                                     std::get<5>( xTup ),
                                     std::get<0>( xTup ),
                                     std::get<6>( xTup ),
                                     std::get<9>( xTup ),
                                     pPack,
                                     pMMIndex,
                                     pRead,
                                     std::get<8>( xTup ),
                                     rParameters.getSelected( )->xRectangularSoc->get( ),
                                     std::get<10>( xTup ) );
                        vColIds.push_back( uiCategoryCounter + ( vEndColumn.size( ) - 1 ) / 2 );

                        addRectangle( vRectangles, xHelperRet, uiCategoryCounter, iReadId, vEndColumn.size( ), false );
                        addRectangle( vRectangles, xReseedOutExtraInfo, uiCategoryCounter, iReadId, vEndColumn.size( ),
                                      true );

                        uiCategoryCounter += vEndColumn.size( ) + 2;
                        vReadsNCols.emplace_back( vColIds.back( ), iReadId );
                        vReads.push_back( pRead );
                    } // else
                } // for
            },
            uiI ) );


    // wait for threads to finish at most iMaxTime seconds, then stop all work and give up
    for( auto& xFuture : vFutures )
    {
        // if maxtime is 0 wait for howeverlong the complete computation takes
        if( iMaxTime > 0 )
        {
            // get status until it is either timeout or ready
            std::future_status xStatus = std::future_status::deferred;
            while( !bStop && xStatus == std::future_status::deferred )
                xStatus = xFuture.wait_until( xEndTime );

            if( xStatus == std::future_status::timeout )
            {
                std::lock_guard<std::mutex> xGuard( xLock );
                bStop = true;
            } // if
        } // if

        xFuture.get( );
    } // for

    if( bStop )
    {
        vRet.clear( );
        vRectangles.clear( );
        vReadsNCols.clear( );
        vReads.clear( );
        vColIds.clear( );
        vAllColIds.clear( );
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
            addSeed( std::get<1>( xTup ), vRet, vEndColumn, vAllColIds, uiCategoryCounter, std::get<3>( xTup ),
                     std::get<4>( xTup ), std::get<2>( xTup ), std::get<5>( xTup ), std::get<0>( xTup ),
                     std::get<6>( xTup ), std::get<9>( xTup ), pPack, pMMIndex, std::get<7>( xTup ),
                     std::get<8>( xTup ), rParameters.getSelected( )->xRectangularSoc->get( ), std::get<10>( xTup ) );
        uiCategoryCounter += vEndColumn.size( );
    } // if

    return std::make_shared<ReadInfo>( ReadInfo{ vRet, vRectangles, vReadsNCols, vReads, vColIds, vAllColIds } );
} // method


#include "ms/container/sv_db/py_db_conf.h"
#include <pybind11/stl.h>
void exportRendererSpeedUp( libMS::SubmoduleOrganizer& xOrganizer )
{
    py::class_<HashCounters, std::shared_ptr<HashCounters>>( xOrganizer.util( ), "HashCounters" ).def( py::init<>( ) );

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
        .def_readwrite( "bOverlapping", &SeedInfo::bOverlapping )
        .def_readwrite( "bInSocReseed", &SeedInfo::bInSocReseed )
        .def_readwrite( "xX", &SeedInfo::xX )
        .def_readwrite( "xY", &SeedInfo::xY )
        .def_readwrite( "uiMinFilterCount", &SeedInfo::uiMinFilterCount )
        .def_readwrite( "uiMaxFilterCount", &SeedInfo::uiMaxFilterCount )
        .def_readwrite( "soc_nt", &SeedInfo::uiSoCNt )
        .def_readwrite( "soc_id", &SeedInfo::uiSocId )
        .def_readwrite( "uiCategory", &SeedInfo::uiCategory );

    py::class_<RectangleInfo>( xOrganizer._util( ), "RectangleInfo" )
        .def_readwrite( "vRectangles", &RectangleInfo::vRectangles )
        .def_readwrite( "vRectangleLayers", &RectangleInfo::vRectangleLayers )
        .def_readwrite( "vRectangleFillPercentage", &RectangleInfo::vRectangleFillPercentage )
        .def_readwrite( "vRectangleReferenceAmbiguity", &RectangleInfo::vRectangleReferenceAmbiguity )
        .def_readwrite( "vRectangleKMerSize", &RectangleInfo::vRectangleKMerSize )
        .def_readwrite( "vRectangleUsedDp", &RectangleInfo::vRectangleUsedDp )
        .def_readwrite( "uiCategory", &RectangleInfo::uiCategory )
        .def_readwrite( "iReadId", &RectangleInfo::iReadId )
        .def_readwrite( "bInSoCReseeding", &RectangleInfo::bInSoCReseeding )
        .def_readwrite( "uiEndColumnSize", &RectangleInfo::uiEndColumnSize );

    py::class_<ReadInfo, std::shared_ptr<ReadInfo>>( xOrganizer._util( ), "ReadInfo" )
        .def_readwrite( "vRet", &ReadInfo::vRet )
        .def_readwrite( "vRectangles", &ReadInfo::vRectangles )
        .def_readwrite( "vReadsNCols", &ReadInfo::vReadsNCols )
        .def_readwrite( "vReads", &ReadInfo::vReads )
        .def_readwrite( "vColIds", &ReadInfo::vColIds )
        .def_readwrite( "vAllColIds", &ReadInfo::vAllColIds );

    xOrganizer.util( ).def( "seedDisplaysForReadIds", &seedDisplaysForReadIds<DBCon> );
} // function
#endif