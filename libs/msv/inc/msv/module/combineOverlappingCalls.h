#pragma once
#include "ms/container/sv_db/pool_container.h"
#include "msv/container/sv_db/query_objects/callInserter.h"

namespace libMSV
{

#define LOG false

/** @brief helper function for combine overlapping calls
 * @details
 * Fetches a call and its supporting jumps.
 * Also computes the average insert size of the supporting jumps.
 */
std::pair<SvCall, size_t>
fetchCall( SQLQuery<DBCon, int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool,
                    std::shared_ptr<CompressedNucSeq>, uint32_t>& xGetCall,
           SQLQuery<DBCon, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, int64_t>& xGetSupportingJumps,
           PriKeyDefaultType iId )
{
    metaMeasureAndLogDuration<LOG>( "xGetCall", [&]( ) { xGetCall.execAndFetch( iId ); } );
    auto xCallTup = xGetCall.get( );
    auto xRet = std::make_pair<SvCall, size_t>( SvCall( std::get<1>( xCallTup ), // uiFromStart
                                                        std::get<2>( xCallTup ), // uiToStart
                                                        std::get<3>( xCallTup ), // uiFromSize
                                                        std::get<4>( xCallTup ), // uiToSize
                                                        std::get<5>( xCallTup ), // from_forward
                                                        std::get<6>( xCallTup ), // to_forward
                                                        std::get<8>( xCallTup ) // supporting_reads
                                                        ),
                                                0 );
    xRet.first.pInsertedSequence = std::get<7>( xCallTup ) == nullptr ? nullptr : std::get<7>( xCallTup )->pUncomNucSeq;
    xRet.first.iId = std::get<0>( xCallTup );

    metaMeasureAndLogDuration<LOG>( "xGetSupprotingJumps",
                                    [&]( ) { xGetSupportingJumps.exec( std::get<0>( xCallTup ) ); } );
    while( xGetSupportingJumps.next( ) )
    {
        auto xSupportTup = xGetSupportingJumps.get( );
        xRet.first.vSupportingJumpIds.push_back( std::get<7>( xSupportTup ) );
        auto pNewJump = std::make_shared<SvJump>( std::get<0>( xSupportTup ), std::get<1>( xSupportTup ),
                                                  std::get<2>( xSupportTup ), std::get<3>( xSupportTup ),
                                                  std::get<4>( xSupportTup ), std::get<5>( xSupportTup ),
                                                  std::get<6>( xSupportTup ), std::get<7>( xSupportTup ) );
        xRet.first.vSupportingJumps.push_back( pNewJump );
        xRet.first.addJumpToEstimateClusterSize( pNewJump );
        xRet.second += xRet.first.vSupportingJumps.back( )->query_distance( );
    } // while
    xRet.second /= xRet.first.vSupportingJumps.size( );
    return xRet;
} // function

/** @brief combine overlapping calls
 * @details
 * Combines all calls that rectangles are overlapping.
 * Multithreads the initial detection of overlapping calls.
 * The actual combination of overlapping calls is single threaded though.
 */
template <typename DBCon>
size_t combineOverlappingCalls( const ParameterSetManager& rParameters,
                                std::shared_ptr<libMS::PoolContainer<DBCon>> pConPool, int64_t iSvCallerId )
{
    size_t uiNumCombines = 0;
    pConPool->xPool.run( [&]( std::shared_ptr<DBCon> pOuterConnection ) {
        SQLStatement<DBCon>( pOuterConnection, "DROP TABLE IF EXISTS call_overlap_table " ).exec( );
        SQLTableWithAutoPriKey<DBCon, PriKeyDefaultType, PriKeyDefaultType> xOverlapTable(
            pOuterConnection,
            json{{TABLE_NAME, "call_overlap_table"},
                 // {CPP_EXTRA, "DROP ON DESTRUCTION"},
                 {TABLE_COLUMNS,
                  {
                      // cannot have references here otherwise msql locks the table if a query is run over sv_call_table
                      {{COLUMN_NAME, "from_call_id"}, /*{REFERENCES, "sv_call_table(id)"},*/ {CONSTRAINTS, "NOT NULL"}},
                      {{COLUMN_NAME, "to_call_id"}, /*{REFERENCES, "sv_call_table(id)"},*/ {CONSTRAINTS, "NOT NULL"}} //
                  }}} );

        metaMeasureAndLogDuration<LOG>( "xInsertIntoOverlapTable", [&]( ) {
            // fill call_overlap_table
            size_t uiNumTasks = rParameters.getNumThreads( ) * 10; // @todo make a parameter out of this
            std::vector<std::future<void>> vFutures;
            vFutures.reserve( uiNumTasks );
            for( size_t uiI = 0; uiI < uiNumTasks; uiI++ )
                vFutures.push_back( pConPool->xPool.enqueue(
                    [&]( std::shared_ptr<DBCon> pConnection, size_t uiI ) {
                        // @todo this triggers multiple full table scans...
                        SQLQuery<DBCon, PriKeyDefaultType, WKBUint64Rectangle, bool, bool> xFetchCalls(
                            pConnection,
                            "SELECT id, ST_AsBinary(rectangle), from_forward, to_forward "
                            "FROM sv_call_table "
                            "WHERE sv_caller_run_id = ? "
                            "AND id % ? = ? ",
                            "combineOverlappingCalls::xFetchCalls" );
                        // split into two queries, so that mysql uses the spatial index for MBRIntersects
                        // stupid query optimizer...
                        SQLStatement<DBCon> xInsertIntoOverlapTable(
                            pConnection,
                            "INSERT INTO call_overlap_table (from_call_id, to_call_id) "
                            "SELECT ?, id "
                            "FROM sv_call_table "
                            "WHERE sv_caller_run_id = ? "
                            "AND id > ? "
                            "AND " ST_INTERSCTS "(rectangle, ST_PolyFromWKB(?, 0)) "
                            "AND from_forward = ? "
                            "AND to_forward = ? ",
                            "combineOverlappingCalls::xInsertIntoOverlapTable" );
                        xFetchCalls.execAndForAll(
                            [&]( PriKeyDefaultType iOuterId, WKBUint64Rectangle& xRect, bool bFromForward,
                                 bool bToForward ) {
                                xInsertIntoOverlapTable.exec( iOuterId, iSvCallerId, iOuterId, xRect, bFromForward,
                                                              bToForward );
                            }, // lambda
                            iSvCallerId, uiNumTasks, uiI );
                    },
                    uiI ) );
            for( size_t uiI = 0; uiI < uiNumTasks; uiI++ )
                vFutures[ uiI ].get( );
        } );

        // create index on overlap table  (mysql cannot deal with empty tables apparently)
        metaMeasureAndLogDuration<LOG>( "xOverlapTable.addIndex", [&]( ) {
            if( SQLQuery<DBCon, uint64_t>( pOuterConnection, "SELECT COUNT(*) FROM call_overlap_table" ).scalar( ) > 0 )
                xOverlapTable.addIndex( json{{INDEX_NAME, "from_to"}, {INDEX_COLUMNS, "from_call_id, to_call_id"}} );
        } );

        // overlap the calls in call_overlap_table
        SQLQuery<DBCon, PriKeyDefaultType> xGetLowestOverlappingCall(
            pOuterConnection,
            "SELECT from_call_id "
            "FROM call_overlap_table "
            "ORDER BY from_call_id ASC "
            "LIMIT 1 ",
            json{}, "combineOverlappingCalls::xGetLowestOverlappingCall" );
        SQLQuery<DBCon, PriKeyDefaultType> xGetOverlappingCalls( pOuterConnection,
                                                                 "SELECT to_call_id "
                                                                 "FROM call_overlap_table "
                                                                 "WHERE from_call_id = ? ",
                                                                 json{},
                                                                 "combineOverlappingCalls::xGetOverlappingCalls" );

        SQLQuery<DBCon, int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, std::shared_ptr<CompressedNucSeq>,
                 uint32_t>
            xGetCall( pOuterConnection,
                      "SELECT id, from_pos, to_pos, from_size, to_size, from_forward, to_forward, inserted_sequence, "
                      "supporting_reads "
                      "FROM sv_call_table "
                      "WHERE id = ? ",
                      json{}, "combineOverlappingCalls::xGetCall" );
        SQLQuery<DBCon, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, int64_t> xGetSupportingJumps(
            pOuterConnection,
            "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, was_mirrored, "
            "sv_jump_table.id "
            "FROM sv_call_support_table "
            "JOIN sv_jump_table ON sv_call_support_table.jump_id = sv_jump_table.id "
            "WHERE sv_call_support_table.call_id = ? ",
            json{}, "combineOverlappingCalls::xGetSupprotingJumps" );
        SQLTable<DBCon, PriKeyDefaultType> xToDeleteTable( pOuterConnection,
                                                           json{{TABLE_NAME, "to_delete_table"},
                                                                {CPP_EXTRA, "DROP ON DESTRUCTION"},
                                                                {TABLE_COLUMNS,
                                                                 {
                                                                     {{COLUMN_NAME, "call_id"},
                                                                      {REFERENCES, "sv_call_table(id)"},
                                                                      {CONSTRAINTS, "NOT NULL PRIMARY KEY"}}
                                                                     //
                                                                 }}} );
        SQLStatement<DBCon> xDoneWithCall( pOuterConnection,
                                           "DELETE FROM call_overlap_table "
                                           "WHERE to_call_id = ? ",
                                           "combineOverlappingCalls::xDoneWithCall" );
        metaMeasureAndLogDuration<LOG>( "xOverlap", [&]( ) {
            while( xGetLowestOverlappingCall.execAndFetch( ) )
            {
                PriKeyDefaultType iLowestCallId = std::get<0>( xGetLowestOverlappingCall.get( ) );
                std::set<PriKeyDefaultType> xIdsToCombineWith;

                // fill xIdsToCombineWith
                std::queue<PriKeyDefaultType> xIdsToCheck;
                xIdsToCheck.push( iLowestCallId );
                while( !xIdsToCheck.empty( ) )
                {
                    xGetOverlappingCalls.execAndForAll(
                        [&]( PriKeyDefaultType iId ) {
                            if( xIdsToCombineWith.find( iId ) == xIdsToCombineWith.end( ) )
                            {
                                xIdsToCombineWith.insert( iId );
                                xIdsToCheck.push( iId );
                            } // if
                        },
                        xIdsToCheck.front( ) );
                    xIdsToCheck.pop( );
                } // while

                auto xPrim = fetchCall( xGetCall, xGetSupportingJumps, iLowestCallId );

                // combine calls
                for( auto iId : xIdsToCombineWith )
                {
                    auto xToCombineWith = fetchCall( xGetCall, xGetSupportingJumps, iId );

                    nucSeqIndex uiMaxInsertSizeDiff = 150;
                    // make sure that the average insert size of the clusters is no more different than
                    // uiMaxInsertSizeDiff this must be done since there could be two different sequences
                    // for the same edge, in that case we should evaluate those sequences independently.
                    if( xPrim.second + uiMaxInsertSizeDiff >= xToCombineWith.second &&
                        xToCombineWith.second + uiMaxInsertSizeDiff >= xPrim.second )
                    {
                        xPrim.first.join( xToCombineWith.first );
                        uiNumCombines++;
                        xToDeleteTable.insert( xToCombineWith.second );
                    } // if
                    xDoneWithCall.exec( xToCombineWith.second );
                } // for
            } // while
        } );

        // delete calls that have been overlapped (mysql cannot deal with empty tables apparently)
        if( SQLQuery<DBCon, uint64_t>( pOuterConnection, "SELECT COUNT(*) FROM to_delete_table" ).scalar( ) > 0 )
        {
            SQLStatement<DBCon> xDeleteOverlapped( pOuterConnection,
                                                   "DELETE FROM sv_call_table "
                                                   "WHERE id IN ( "
                                                   "    SELECT call_id "
                                                   "    FROM to_delete_table "
                                                   ") ",
                                                   "combineOverlappingCalls::xDeleteOverlapped" );

            metaMeasureAndLogDuration<LOG>( "xDeleteOverlapped", [&]( ) { xDeleteOverlapped.exec( ); } );
        } // if
    } ); // xPool.run call

    return uiNumCombines;

} // function
}; // namespace libMSV
