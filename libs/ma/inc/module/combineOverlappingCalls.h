#pragma once
#include "container/sv_db/query_objects/callInserter.h"

namespace libMA
{
// @todo this can be done via a linesweep...
template <typename DBCon>
size_t combineOverlappingCalls( const ParameterSetManager& rParameters, std::shared_ptr<SV_Schema<DBCon>> pDb,
                                int64_t iSvCallerId )
{
    SQLQuery<DBCon, int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, NucSeqSql, uint32_t> xQuery(
        pDb->pDatabase,
        "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, inserted_sequence, supporting_reads "
        "FROM sv_call_table "
        "WHERE sv_caller_run_id = ? "
        "ORDER BY id " );

    SQLQuery<DBCon, int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, NucSeqSql, uint32_t> xQuery2(
        pDb->pDatabase,
        "SELECT sv_call_r_tree.id, from_pos, to_pos, from_size, to_size, switch_strand, "
        "       inserted_sequence, supporting_reads "
        "FROM sv_call_table, sv_call_r_tree "
        "WHERE sv_call_table.id = sv_call_r_tree.id "
        "AND sv_call_r_tree.run_id_a >= ? " // dim 1
        "AND sv_call_r_tree.run_id_b <= ? " // dim 1
        "AND sv_call_r_tree.maxX >= ? " // dim 2
        "AND sv_call_r_tree.minX <= ? " // dim 2
        "AND sv_call_r_tree.maxY >= ? " // dim 3
        "AND sv_call_r_tree.minY <= ? " // dim 3
        "AND sv_call_r_tree.id > ? "
        "AND switch_strand = ? " );

    SQLQuery<DBCon, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, int64_t> xQuerySupport(
        pDb->pDatabase,
        "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
        "sv_jump_table.id "
        "FROM sv_call_support_table "
        "JOIN sv_jump_table ON sv_call_support_table.jump_id = sv_jump_table.id "
        "WHERE sv_call_support_table.call_id = ? " );

    SvCallInserter<DBCon> xInserter( pDb, iSvCallerId ); // also triggers transaction
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
            // DEL: auto xIt2 = xQuery2.vExecuteAndReturnIterator( iSvCallerId, iSvCallerId,
            // DEL:                                                std::get<1>( xTup ), // from_pos
            // DEL:                                                // from_pos + from_size
            // DEL:                                                std::get<1>( xTup ) + std::get<3>( xTup ),
            // DEL:                                                std::get<2>( xTup ), // to_pos
            // DEL:                                                // to_pos + to_size
            // DEL:                                                std::get<2>( xTup ) + std::get<4>( xTup ),
            // DEL:                                                std::get<0>( xTup ), // id
            // DEL:                                                std::get<5>( xTup ) // bSwitchStrand
            // DEL: );
            xQuery2.execAndFetch( iSvCallerId, iSvCallerId,
                                  std::get<1>( xTup ), // from_pos
                                  std::get<1>( xTup ) + std::get<3>( xTup ), // from_pos + from_size
                                  std::get<2>( xTup ), // to_pos
                                  std::get<2>( xTup ) + std::get<4>( xTup ), // to_pos + to_size
                                  std::get<0>( xTup ), // id
                                  std::get<5>( xTup ) ); // bSwitchStrand)
            // DEL: if( !xIt2.eof( ) )
            if( !xQuery2.eof( ) ) // FIXME: avoid using eof - work with next() and get() merely
            {

                xPrim.pInsertedSequence = std::get<6>( xTup ).pNucSeq;
                xPrim.iId = std::get<0>( xTup );
                // DEL: auto xSupportIterator( xQuerySupport.vExecuteAndReturnIterator( std::get<0>( xTup ) ) );
                xQuerySupport.execAndFetch( std::get<0>( xTup ) );
                nucSeqIndex uiPrimInsertSizeAvg = 0;
                // DEL: while( !xSupportIterator.eof( ) )
                while( !xQuerySupport.eof( ) )
                {
                    // DEL: auto xTup = xSupportIterator.get( );
                    auto xTup = xQuerySupport.get( );
                    xPrim.vSupportingJumpIds.push_back( std::get<7>( xTup ) );
                    auto pNewJump =
                        std::make_shared<SvJump>( rParameters.getSelected( ), std::get<0>( xTup ), std::get<1>( xTup ),
                                                  std::get<2>( xTup ), std::get<3>( xTup ), std::get<4>( xTup ),
                                                  std::get<5>( xTup ), std::get<6>( xTup ), std::get<7>( xTup ) );
                    xPrim.vSupportingJumps.push_back( pNewJump );
                    xPrim.addJumpToEstimateClusterSize( pNewJump );
                    uiPrimInsertSizeAvg += xPrim.vSupportingJumps.back( )->query_distance( );
                    // DEL: xSupportIterator.next( );
                    xQuerySupport.next( );
                } // while
                uiPrimInsertSizeAvg /= xPrim.vSupportingJumps.size( );


                // DEL: while( !xIt2.eof( ) )
                while( !xQuery2.eof( ) )
                {
                    // DEL: auto xTup2 = xIt2.get( );
                    auto xTup2 = xQuery2.get( );
                    // DEL: xIt2.next( );
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
                    // DEL: auto xSupportIterator( xQuerySupport.vExecuteAndReturnIterator( std::get<0>( xTup2 ) ) );
                    xQuerySupport.execAndFetch( std::get<0>( xTup2 ) );
                    nucSeqIndex uiSecInsertSizeAvg = 0;
                    // DEL: while( !xSupportIterator.eof( ) )
                    while( !xQuerySupport.eof( ) )
                    {
                        // DEL: auto xTup = xSupportIterator.get( );
                        auto xTup = xQuerySupport.get( );
                        xSec.vSupportingJumpIds.push_back( std::get<7>( xTup ) );
                        auto pNewJump = std::make_shared<SvJump>(
                            rParameters.getSelected( ), std::get<0>( xTup ), std::get<1>( xTup ), std::get<2>( xTup ),
                            std::get<3>( xTup ), std::get<4>( xTup ), std::get<5>( xTup ), std::get<6>( xTup ),
                            std::get<7>( xTup ) );
                        xSec.vSupportingJumps.push_back( pNewJump );
                        xSec.addJumpToEstimateClusterSize( pNewJump );
                        uiSecInsertSizeAvg += xSec.vSupportingJumps.back( )->query_distance( );
                        // DEL:xSupportIterator.next( );
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
                pDb->pSvCallSupportTable->deleteCall( iId );
                pDb->pSvCallTable->deleteCall( iId );
                sDeleted.insert( iId );
                uiRet++;
            } // for
            xPrim.reEstimateClusterSize( );
            xInserter.updateCall( xPrim );
        } // if

    } // for

    return uiRet;
    // end of scope for transaction context (via xInserter)

} // function
}; // namespace libMA
