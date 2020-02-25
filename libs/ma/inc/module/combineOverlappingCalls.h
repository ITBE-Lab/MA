#pragma once
#include "container/sv_db/query_objects/callInserter.h"

namespace libMA
{
// @todo this can be done via a linesweep...
template <typename DBCon> size_t combineOverlappingCalls( std::shared_ptr<DBCon> pConnection, int64_t iSvCallerId )
{
    ExplainedSQLQuery<DBCon, int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, NucSeqSql, uint32_t> xQuery(
        pConnection,
        "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, inserted_sequence, supporting_reads "
        "FROM sv_call_table "
        "WHERE sv_caller_run_id = ? "
        "ORDER BY id ",
        json{}, "combineOverlappingCalls::xQuery" );

    ExplainedSQLQuery<DBCon, int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, NucSeqSql, uint32_t> xQuery2(
        pConnection,
        "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, "
        "       inserted_sequence, supporting_reads "
        "FROM sv_call_table "
        "WHERE sv_caller_run_id = ? "
        "AND MBRIntersects(rectangle, ST_PolyFromWKB(?, 0)) "
        "AND id > ? "
        "AND switch_strand = ? ",
        json{}, "combineOverlappingCalls::xQuery2" );

    ExplainedSQLQuery<DBCon, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, int64_t> xQuerySupport(
        pConnection,
        "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
        "sv_jump_table.id "
        "FROM sv_call_support_table "
        "JOIN sv_jump_table ON sv_call_support_table.jump_id = sv_jump_table.id "
        "WHERE sv_call_support_table.call_id = ? ",
        json{}, "combineOverlappingCalls::xQuerySupport" );

    SvCallTable<DBCon> xSvTable( pConnection );
    SvCallSupportTable<DBCon> xSvCallSupportTable( pConnection );

    auto vQuery1Res = xQuery.executeAndStoreAllInVector( iSvCallerId );

    size_t uiRet = 0;
    // iterate over all calls

    std::set<int64_t> sDeleted;
    for( auto& xTup : vQuery1Res )
    {
        if( sDeleted.find( std::get<0>( xTup ) ) != sDeleted.end( ) )
            continue;

        SvCall xPrim( std::get<1>( xTup ), // uiFromStart
                      std::get<2>( xTup ), // uiToStart
                      std::get<3>( xTup ), // uiFromSize
                      std::get<4>( xTup ), // uiToSize
                      std::get<5>( xTup ), // bSwitchStrand
                      std::get<7>( xTup ) // supporting_reads
        );

        std::vector<int64_t> vToDel;
        {
            // get all overlapping calls

            auto xWkb = WKBUint64Rectangle( geom::Rectangle<nucSeqIndex>( std::get<1>( xTup ), std::get<2>( xTup ),
                                                                          std::get<3>( xTup ), std::get<4>( xTup ) ) );
            xQuery2.execAndFetch( iSvCallerId,
                                  xWkb, // rectangle
                                  std::get<0>( xTup ), // id
                                  std::get<5>( xTup ) ); // bSwitchStrand)
            if( !xQuery2.eof( ) ) // FIXME: avoid using eof - work with next() and get() merely
            {

                xPrim.pInsertedSequence = std::get<6>( xTup ).pNucSeq;
                xPrim.iId = std::get<0>( xTup );
                xQuerySupport.execAndFetch( std::get<0>( xTup ) );
                nucSeqIndex uiPrimInsertSizeAvg = 0;
                while( !xQuerySupport.eof( ) )
                {
                    auto xTup = xQuerySupport.get( );
                    xPrim.vSupportingJumpIds.push_back( std::get<7>( xTup ) );
                    auto pNewJump = std::make_shared<SvJump>(
                        std::get<0>( xTup ), std::get<1>( xTup ), std::get<2>( xTup ), std::get<3>( xTup ),
                        std::get<4>( xTup ), std::get<5>( xTup ), std::get<6>( xTup ), std::get<7>( xTup ) );
                    xPrim.vSupportingJumps.push_back( pNewJump );
                    xPrim.addJumpToEstimateClusterSize( pNewJump );
                    uiPrimInsertSizeAvg += xPrim.vSupportingJumps.back( )->query_distance( );
                    // DEL: xSupportIterator.next( );
                    xQuerySupport.next( );
                } // while
                uiPrimInsertSizeAvg /= xPrim.vSupportingJumps.size( );


                while( !xQuery2.eof( ) )
                {
                    auto xTup2 = xQuery2.get( );
                    xQuery2.next( );

                    if( sDeleted.find( std::get<0>( xTup2 ) ) != sDeleted.end( ) )
                        continue;

                    SvCall xSec( std::get<1>( xTup2 ), // uiFromStart
                                 std::get<2>( xTup2 ), // uiToStart
                                 std::get<3>( xTup2 ), // uiFromSize
                                 std::get<4>( xTup2 ), // uiToSize
                                 std::get<5>( xTup2 ), // bSwitchStrand
                                 std::get<7>( xTup2 ) // supporting_reads
                    );
                    xSec.pInsertedSequence = std::get<6>( xTup2 ).pNucSeq;
                    xSec.iId = std::get<0>( xTup2 );
                    xQuerySupport.execAndFetch( std::get<0>( xTup2 ) );
                    nucSeqIndex uiSecInsertSizeAvg = 0;
                    while( !xQuerySupport.eof( ) )
                    {
                        auto xTup = xQuerySupport.get( );
                        xSec.vSupportingJumpIds.push_back( std::get<7>( xTup ) );
                        auto pNewJump = std::make_shared<SvJump>(
                            std::get<0>( xTup ), std::get<1>( xTup ), std::get<2>( xTup ), std::get<3>( xTup ),
                            std::get<4>( xTup ), std::get<5>( xTup ), std::get<6>( xTup ), std::get<7>( xTup ) );
                        xSec.vSupportingJumps.push_back( pNewJump );
                        xSec.addJumpToEstimateClusterSize( pNewJump );
                        uiSecInsertSizeAvg += xSec.vSupportingJumps.back( )->query_distance( );
                        xQuerySupport.next( );
                    } // while
                    uiSecInsertSizeAvg /= xSec.vSupportingJumps.size( );

                    nucSeqIndex uiMaxInsertSizeDiff = 150;
                    // make sure that the average insert size of the clusters is no more different than
                    // uiMaxInsertSizeDiff this must be done since there could be two different sequences for the same
                    // edge, in that case we should evaluate those sequences independently.
                    if( uiPrimInsertSizeAvg + uiMaxInsertSizeDiff >= uiSecInsertSizeAvg &&
                        uiSecInsertSizeAvg + uiMaxInsertSizeDiff >= uiPrimInsertSizeAvg )
                    {
                        xPrim.join( xSec );

                        vToDel.push_back( xSec.iId );
                    } // if
                } // while
            } // if
        } // scope for xIt2

        // only replace the primary cluster if we actually made changes
        if( vToDel.size( ) > 0 )
        {
            for( int64_t iId : vToDel )
            {
                xSvCallSupportTable.deleteCall( iId );
                xSvTable.deleteCall( iId );
                sDeleted.insert( iId );
                uiRet++;
            } // for
            xPrim.reEstimateClusterSize( );
            xSvTable.updateCall( iSvCallerId, xPrim );
        } // if

    } // for

    return uiRet;
    // end of scope for transaction context (via xInserter)

} // function
}; // namespace libMA
