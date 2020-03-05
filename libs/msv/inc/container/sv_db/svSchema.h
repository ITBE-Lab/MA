/**
 * @file svDb.h
 * @details
 * The database interface for the structural variant caller
 */
#pragma once

#include "container/pack.h"
#include "container/sv_db/pool_container.h"
#include "module/module.h"

namespace libMA
{

template <typename DBCon> class CountCallsQuery : public SQLQuery<DBCon, uint32_t>
{
  public:
    CountCallsQuery( std::shared_ptr<DBCon> pConn )
        : SQLQuery<DBCon, uint32_t>( pConn, "SELECT COUNT(*) "
                                            "FROM sv_call_table "
                                            "WHERE sv_caller_run_id = ? "
                                            "AND MBRIntersects(rectangle, ST_PolyFromWKB(?, 0)) "
                                            "AND score >= ? "
                                            "LIMIT ? " )
    {} // constructor
}; // class

template <typename DBCon>
uint32_t getCallOverviewAreaHelper( CountCallsQuery<DBCon>& xQuery, std::shared_ptr<Pack> pPack, int64_t iRunId,
                                    double dMinScore, int64_t iX, int64_t iY, uint64_t uiW, uint64_t uiH,
                                    uint32_t uiLimit )
{
    uint32_t uiX = 0;
    if( iX > 0 )
        uiX = (uint32_t)iX;
    uint32_t uiY = 0;
    if( iY > 0 )
        uiY = (uint32_t)iY;
    if( uiX + uiW > pPack->uiUnpackedSizeForwardStrand )
        uiW = pPack->uiUnpackedSizeForwardStrand - uiX;
    if( uiY + uiH > pPack->uiUnpackedSizeForwardStrand )
        uiH = pPack->uiUnpackedSizeForwardStrand - uiY;

    auto xWkb = WKBUint64Rectangle( geom::Rectangle<nucSeqIndex>( uiX, uiY, uiW, uiH ) );
    return xQuery.scalar( iRunId, xWkb, dMinScore, uiLimit );
} // function

template <typename DBCon>
uint32_t getCallOverviewArea( std::shared_ptr<DBCon> pConnection, std::shared_ptr<Pack> pPack, int64_t iRunId,
                              double dMinScore, int64_t iX, int64_t iY, uint64_t uiW, uint64_t uiH, uint32_t uiLimit )
{
    CountCallsQuery<DBCon> xQuery( pConnection );
    return getCallOverviewAreaHelper( xQuery, pPack, iRunId, dMinScore, iX, iY, uiW, uiH, uiLimit );
} // function


struct rect
{
    uint32_t x, y, w, h, c, i, j;
    rect( uint32_t x, uint32_t y, uint32_t w, uint32_t h, uint32_t c, uint32_t i, uint32_t j )
        : x( x ), y( y ), w( w ), h( h ), c( c ), i( i ), j( j )
    {} // constructor
}; // struct

template <typename DBCon>
std::vector<rect> getCallOverview( std::shared_ptr<PoolContainer<DBCon>> pConPool, std::shared_ptr<Pack> pPack,
                                   int64_t iRunId, double dMinScore, int64_t iX, int64_t iY, uint64_t uiW, uint64_t uiH,
                                   uint64_t uiMaxW, uint64_t uiMaxH, uint32_t uiGiveUpFactor )
{
    uint32_t uiX = 0;
    if( iX > 0 )
        uiX = (uint32_t)iX;
    uint32_t uiY = 0;
    if( iY > 0 )
        uiY = (uint32_t)iY;
    if( uiX + uiW > pPack->uiUnpackedSizeForwardStrand )
        uiW = pPack->uiUnpackedSizeForwardStrand - uiX;
    if( uiY + uiH > pPack->uiUnpackedSizeForwardStrand )
        uiH = pPack->uiUnpackedSizeForwardStrand - uiY;
    size_t uiStartContigX = pPack->uiSequenceIdForPosition( uiX );
    size_t uiEndContigX = pPack->uiSequenceIdForPosition( uiX + uiW );
    size_t uiStartContigY = pPack->uiSequenceIdForPosition( uiY );
    size_t uiEndContigY = pPack->uiSequenceIdForPosition( uiY + uiH );
    std::vector<rect> vRet;

    std::vector<std::unique_ptr<CountCallsQuery<DBCon>>> vQueries;
    // vector size cannot be adjusted in parallel; so this is necessary
    vQueries.resize( pConPool->xPool.uiPoolSize );

    std::vector<std::future<std::vector<rect>>> vFutures;

    for( size_t uiContigX = uiStartContigX; uiContigX <= uiEndContigX; uiContigX++ )
        for( size_t uiContigY = uiStartContigY; uiContigY <= uiEndContigY; uiContigY++ )
        {
            auto uiStartX = std::max( uiX, (uint32_t)pPack->startOfSequenceWithId( uiContigX ) );
            auto uiEndX = std::min( uiX + uiW, (uint64_t)pPack->endOfSequenceWithId( uiContigX ) );
            auto uiStartY = std::max( uiY, (uint32_t)pPack->startOfSequenceWithId( uiContigY ) );
            auto uiEndY = std::min( uiY + uiH, (uint64_t)pPack->endOfSequenceWithId( uiContigY ) );
            uint32_t uiNumW = ( uint32_t )( uiEndX - uiStartX ) / (uint32_t)uiMaxW + (uint32_t)1;
            uint32_t uiNumH = ( uint32_t )( uiEndY - uiStartY ) / (uint32_t)uiMaxH + (uint32_t)1;
            double dW = ( (double)( uiEndX - uiStartX ) ) / (double)uiNumW;
            double dH = ( (double)( uiEndY - uiStartY ) ) / (double)uiNumH;
            if( dW * uiGiveUpFactor < uiW )
                continue;
            if( dH * uiGiveUpFactor < uiH )
                continue;
            for( size_t uiI = 0; uiI < uiNumW; uiI++ )
                vFutures.push_back( pConPool->xPool.enqueue(
                    [&]( std::shared_ptr<DBCon> pConnection, size_t uiContigX, size_t uiContigY, size_t uiI, double dW,
                         double dH, uint32_t uiStartX, uint32_t uiStartY, uint32_t uiNumH ) {
                        std::vector<rect> vRet;
                        for( size_t uiJ = 0; uiJ < uiNumH; uiJ++ )
                        {
                            // init connection if necessary
                            if( vQueries[ pConnection->getTaskId( ) ] == nullptr )
                                vQueries[ pConnection->getTaskId( ) ] =
                                    std::make_unique<CountCallsQuery<DBCon>>( pConnection );

                            int64_t uiInnerX = ( int64_t )( uiI * dW + uiStartX );
                            int64_t uiInnerY = ( int64_t )( uiJ * dH + uiStartY );
                            auto c = getCallOverviewAreaHelper(
                                *vQueries[ pConnection->getTaskId( ) ], pPack, iRunId, dMinScore, uiInnerX, uiInnerY,
                                (uint32_t)dW + 1, (uint32_t)dH + 1, std::numeric_limits<uint32_t>::max( ) );
                            if( c > 0 )
                                vRet.emplace_back( (uint32_t)uiInnerX, (uint32_t)uiInnerY, (uint32_t)dW, (uint32_t)dH,
                                                   c, (uint32_t)uiContigX, (uint32_t)uiContigY );
                        } // for
                        return vRet;
                    },
                    uiContigX, uiContigY, uiI, dW, dH, uiStartX, uiStartY, uiNumH ) );
        } // for
    for( auto& xFuture : vFutures )
        for( auto& xRect : xFuture.get( ) )
            vRet.push_back( xRect );
    return vRet;
} // function

} // namespace libMA

#ifdef WITH_PYTHON
void exportSoCDbWriter( py::module& rxPyModuleId );
#endif
