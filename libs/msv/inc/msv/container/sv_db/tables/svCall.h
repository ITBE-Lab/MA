/**
 * @file svCall.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "ma/container/pack.h"
#include "ms/container/sv_db/pool_container.h"
#include "msv/container/svJump.h"
#include "msv/container/sv_db/tables/svCallerRun.h"
#include "sql_api.h"
#include "util/geom.h"
#include "util/system.h"
#include "util/threadPool.h"
#include "wkb_spatial.h"
#include <csignal>
#include <string>
#include <memory>
#include <type_traits>

using namespace libMA;
namespace libMSV
{


template <typename DBCon>
using SvCallTableType = SQLTableWithLibIncrPriKey<DBCon, // DB connector type
                                                  PriKeyDefaultType, // sv_caller_run_id (foreign key)
                                                  uint32_t, // from_pos (geometry)
                                                  uint32_t, // to_pos (geometry)
                                                  uint32_t, // from_size (geometry)
                                                  uint32_t, // to_size (geometry)
                                                  bool, // from_forward
                                                  bool, // to_forward
                                                  std::shared_ptr<CompressedNucSeq>, // inserted_sequence
                                                  uint32_t, // inserted_sequence_size
                                                  uint32_t, // supporting_reads
                                                  uint32_t, // supporting_nt
                                                  uint32_t, // reference_ambiguity
                                                  int64_t, // order_id
                                                  int64_t, // ctg_order_id
                                                  bool, // mirrored
                                                  WKBUint64Rectangle // rectangle (geometry)
                                                  >;

template <typename DBCon> class SvCallTable : public SvCallTableType<DBCon>
{
    std::shared_ptr<DBCon> pConnection;
    SQLQuery<DBCon, uint64_t> xQuerySize;
    SQLQuery<DBCon, uint64_t> xQuerySizeSpecific;
    SQLQuery<DBCon, int64_t> xCallArea;
    SQLQuery<DBCon, double> xMaxScore;
    SQLQuery<DBCon, double> xMinScore;
    SQLQuery<DBCon, double> xMaxAvSuppNt;
    SQLStatement<DBCon> xSetCoverageForCall;
    SQLStatement<DBCon> xVacuumAnalyze;
    SQLStatement<DBCon> xDeleteCall;
    SQLStatement<DBCon> xUpdateCall;
    SQLStatement<DBCon> xFilterCallsWithHighScore;
    SQLStatement<DBCon> xEnableExtension;
    SQLStatement<DBCon> xCopyPath;
    SQLStatement<DBCon> xExtractSmallCalls;

  public:
    // Consider: Place the table on global level
    json jSvCallTableDef( )
    {
        return json{ { TABLE_NAME, "sv_call_table" },
                     { TABLE_COLUMNS,
                       { { { COLUMN_NAME, "sv_caller_run_id" } },
                         { { COLUMN_NAME, "from_pos" } },
                         { { COLUMN_NAME, "to_pos" } },
                         { { COLUMN_NAME, "from_size" } },
                         { { COLUMN_NAME, "to_size" } },
                         { { COLUMN_NAME, "from_forward" } },
                         { { COLUMN_NAME, "to_forward" } },
                         { { COLUMN_NAME, "inserted_sequence" } },
                         { { COLUMN_NAME, "inserted_sequence_size" } },
                         { { COLUMN_NAME, "supporting_reads" } },
                         { { COLUMN_NAME, "supporting_nt" } },
                         { { COLUMN_NAME, "reference_ambiguity" } },
                         { { COLUMN_NAME, "order_id" } },
                         { { COLUMN_NAME, "ctg_order_id" } },
                         { { COLUMN_NAME, "mirrored" } },
                         { { COLUMN_NAME, "rectangle" }, { CONSTRAINTS, "NOT NULL" } } } },
                     { GENERATED_COLUMNS,
                       { { { COLUMN_NAME, "score" },
                           { TYPE, DBCon::TypeTranslator::template getSQLColumnTypeName<double>( ) },
                           { AS, "( supporting_reads * 1.0 ) / reference_ambiguity" } },
                         { { COLUMN_NAME, "flipped_rectangle" },
                           { TYPE, DBCon::TypeTranslator::template getSQLColumnTypeName<WKBUint64Rectangle>( ) },
                           { AS, "ST_FlipCoordinates(rectangle::geometry)" } },
                         { { COLUMN_NAME, "avg_supporting_nt" },
                           { TYPE, DBCon::TypeTranslator::template getSQLColumnTypeName<double>( ) },
                           { AS, "( supporting_nt * 1.0 ) / GREATEST(supporting_reads * 1.0, 1.0)" } } } } };
    }; // namespace libMSV

    SvCallTable( std::shared_ptr<DBCon> pConnection )
        : SvCallTableType<DBCon>( pConnection, // the database where the table resides
                                  jSvCallTableDef( ) ), // table definition
          pConnection( pConnection ),
          xQuerySize( pConnection, "SELECT COUNT(*) FROM sv_call_table" ),
          xQuerySizeSpecific( pConnection, "SELECT COUNT(*) FROM sv_call_table "
                                           "WHERE sv_caller_run_id = ? "
                                           "AND score >= ? " ),
          xCallArea( pConnection,
                     "SELECT SUM( from_size * to_size ) FROM sv_call_table "
                     "WHERE sv_caller_run_id = ? "
                     "AND score >= ? " ),
          xMaxScore( pConnection,
                     "SELECT score "
                     " FROM sv_call_table "
                     "WHERE sv_caller_run_id = ? "
                     "ORDER BY score DESC LIMIT 1 " ),
          xMinScore( pConnection,
                     "SELECT score "
                     " FROM sv_call_table "
                     "WHERE sv_caller_run_id = ? "
                     "ORDER BY score ASC LIMIT 1 " ),
          xMaxAvSuppNt( pConnection,
                        "SELECT avg_supporting_nt "
                        " FROM sv_call_table "
                        "WHERE sv_caller_run_id = ? "
                        "ORDER BY avg_supporting_nt DESC LIMIT 1 " ),
          xSetCoverageForCall( pConnection,
                               "UPDATE sv_call_table "
                               "SET reference_ambiguity = ? "
                               "WHERE id = ? " ),
          xVacuumAnalyze( pConnection, "VACUUM ANALYZE sv_call_table " ),
          xDeleteCall( pConnection,
                       "DELETE FROM sv_call_table "
                       "WHERE id = ? " ),
          xUpdateCall( pConnection,
                       "UPDATE sv_call_table "
                       "SET from_pos = ?, "
                       "    to_pos = ?, "
                       "    from_size = ?, "
                       "    to_size = ?, "
                       "    from_forward = ?, "
                       "    to_forward = ?, "
                       "    inserted_sequence = ?, "
                       "    inserted_sequence_size = ?, "
                       "    supporting_reads = ?, "
                       "    supporting_nt = ?, "
                       "    reference_ambiguity = ?, "
                       "    order_id = ?, "
                       "    ctg_order_id = ?, "
                       "    mirrored = ?, "
                       "    rectangle = ST_GeomFromWKB(?, 0) "
                       "WHERE id = ? " ),
          xFilterCallsWithHighScore( pConnection,
                                     "DELETE FROM sv_call_table "
                                     "WHERE sv_caller_run_id = ? "
                                     // filters out ?% of the calls with the highest scores in sv_call_table
                                     "AND score >= ? " ),
          xEnableExtension( pConnection, "CREATE EXTENSION IF NOT EXISTS btree_gist" ),
          // https://www.postgresql.org/docs/current/sql-update.html
          // Query does not update anything where the FROM & WHERE section returns nothing
          // "query may give unexpected results if" [one row is updated multiple times]
          // However: our query should prevent that as long as the allowed distances are set properly
          xCopyPath( pConnection,
                     std::string( //
                         "UPDATE sv_call_table AS outer_ "
                         "SET order_id=subquery.order_id, "
                         "    inserted_sequence=subquery.inserted_sequence, "
                         "    ctg_order_id=subquery.ctg_order_id, "
                         "    mirrored=subquery.mirrored, "
                         "    inserted_sequence_size=subquery.inserted_sequence_size "
                         "FROM ( "
                         "       SELECT id, order_id, ctg_order_id, inserted_sequence, inserted_sequence_size, "
                         "              rectangle, from_forward, to_forward, mirrored "
                         "       FROM sv_call_table AS inner_ "
                         "       WHERE sv_caller_run_id = ? " ) +
                         //      enforce that subquery has no overlap with it's own run
                         /*   */ selfIntersectionSQL( "inner_", "inner_2" ) +
                         "     ) AS subquery "
                         "WHERE outer_.sv_caller_run_id = ? " +
                         // enforce that outer and inner overlap
                         rectanglesOverlapSQL( "outer_", "subquery" ) +
                         // enforce that outer has no overlap with it's own run
                         selfIntersectionSQL( "outer_", "outer_2" ) ),
          xExtractSmallCalls( pConnection,
                                     "UPDATE sv_call_table "
                                     "SET sv_caller_run_id = ? "
                                     "WHERE sv_caller_run_id = ? "
                                     "AND GREATEST(to_pos - from_pos, inserted_sequence_size) < ? " )
    {} // default constructor


    static std::string rectanglesOverlapSQL_1( std::string sFromTable, std::string sToTable )
    {
        // ST_DWithin returns true if the geometries are within the specified distance of one another
        // makes use of indices
        return std::string( "ST_DWithin(" ) + sToTable + ".rectangle::geometry, " + sFromTable +
               ( sFromTable.empty( ) ? "" : "." ) + "rectangle::geometry, ?) AND " + sToTable +
               ".from_forward = " + sFromTable + ( sFromTable.empty( ) ? "" : "." ) + "from_forward AND " + sToTable +
               ".to_forward = " + sFromTable + ( sFromTable.empty( ) ? "" : "." ) + "to_forward ";
    }
    static std::string rectanglesOverlapSQL_2( std::string sFromTable, std::string sToTable )
    {
        // ST_DWithin returns true if the geometries are within the specified distance of one another
        // makes use of indices
        return std::string( "ST_DWithin(" ) + sToTable +
               ".rectangle::geometry, "
               // this considers reverse complemented calls
               // Note: we need to mirror the rectangle on the matrix diagonal
               //       and invert the strand information
               + sFromTable + ( sFromTable.empty( ) ? "" : "." ) + "flipped_rectangle::geometry, ?) AND " + sToTable +
               ".from_forward != " + sFromTable + ( sFromTable.empty( ) ? "" : "." ) + "to_forward AND " + sToTable +
               ".to_forward != " + sFromTable + ( sFromTable.empty( ) ? "" : "." ) + "from_forward ";
    }
    static std::string rectanglesOverlapSQL( std::string sFromTable, std::string sToTable )
    {
        return std::string( "AND (( " ) + rectanglesOverlapSQL_1( sFromTable, sToTable ) + " ) " + //
               "                 OR ( " + rectanglesOverlapSQL_2( sFromTable, sToTable ) + " )) ";
    }

    static std::string selfIntersectionSQL( std::string sFromTable, std::string sToTable )
    {
        // clang-format off
        std::string sPref = std::string("AND NOT EXISTS( " ) +
            // make sure that inner_table does not overlap with any other call with higher score
            "     SELECT " + sToTable + ".id "
            "     FROM sv_call_table AS " + sToTable +
            "     WHERE " + sToTable + ".id != " + sFromTable + ( sFromTable.empty( ) ? "" : "." ) + "id "
            "     AND " + sToTable + ".score >= " + sFromTable + ( sFromTable.empty( ) ? "" : "." ) + "score "
            "     AND " + sToTable + ".sv_caller_run_id = " + sFromTable + ( sFromTable.empty( ) ? "" : "." ) + 
                                                                                                    "sv_caller_run_id "
            "     AND ";
        // clang-format on
        return sPref + rectanglesOverlapSQL_1( sFromTable, sToTable ) + ") " + sPref +
               rectanglesOverlapSQL_2( sFromTable, sToTable ) + ") ";
    }

    inline void genIndices( int64_t iCallerRunId )
    {
        xEnableExtension.exec( );
        this->addIndex( json{ { INDEX_NAME, "rectangle" },
                              { INDEX_COLUMNS, "sv_caller_run_id, rectangle" },
                              { INDEX_TYPE, "SPATIAL" },
                              { INCLUDE, "id, from_forward, to_forward, rectangle, score" },
                              { INDEX_METHOD, "GIST" } } );
        // see: https://dev.mysql.com/doc/refman/5.7/en/create-table-generated-columns.html
        // and: https://dev.mysql.com/doc/refman/5.7/en/create-table-secondary-indexes.html
        this->addIndex( json{ { INDEX_NAME, "runId_score" }, { INDEX_COLUMNS, "sv_caller_run_id, score" } } );
        this->addIndex( json{ { INDEX_NAME, "from_pos_idx" }, { INDEX_COLUMNS, "from_pos" } } );
        this->addIndex( json{ { INDEX_NAME, "to_pos_idx" }, { INDEX_COLUMNS, "to_pos" } } );
        this->addIndex( json{ { INDEX_NAME, "reconstruction_index" },
                              { INDEX_COLUMNS, "ctg_order_id, sv_caller_run_id, order_id" } } );
    } // method

    inline void dropIndices( int64_t iCallerRunId )
    {
        this->dropIndex( json{ { INDEX_NAME, "rectangle" } } );
        this->dropIndex( json{ { INDEX_NAME, "runId_score" } } );
        this->dropIndex( json{ { INDEX_NAME, "from_pos_idx" } } );
        this->dropIndex( json{ { INDEX_NAME, "to_pos_idx" } } );
        this->dropIndex( json{ { INDEX_NAME, "reconstruction_index" } } );
    } // method

    inline void vacuumAnalyze( )
    {
        xVacuumAnalyze.exec( );
    } // method

    inline uint32_t numCalls( )
    {
        return xQuerySize.scalar( );
    } // method

    inline uint32_t numCalls( int64_t iCallerRunId, double dMinScore )
    {
        return xQuerySizeSpecific.scalar( iCallerRunId, dMinScore );
    } // method
    inline uint32_t numCalls_py( int64_t iCallerRunId, double dMinScore )
    {
        return numCalls( iCallerRunId, dMinScore );
    } // method

    inline void updateCoverage( SvCall& rCall )
    {
        xSetCoverageForCall.exec( (uint32_t)rCall.uiReferenceAmbiguity, rCall.iId );
    } // method

    inline void deleteCall( int64_t iCallId )
    {
        xDeleteCall.exec( iCallId );
    } // method

    inline void deleteCall( SvCall& rCall )
    {
        deleteCall( rCall.iId );
    } // method

    inline void extractSmallCalls(int64_t iCallerRunId, int64_t iMaxSize, std::string sName, std::string sDesc)
    {
        auto pRun = std::make_shared<SvCallerRunTable<DBCon>>( pConnection );
        auto iNewId = pRun->insert( sName, sDesc, pRun->getSvJumpRunId( iCallerRunId) );
        xExtractSmallCalls.exec( iNewId, iCallerRunId, iMaxSize );
        genIndices( iNewId );
    } // method

    inline void copyPath( int64_t iCallerRunIdFrom, int64_t iCallerRunIdTo, int64_t iAllowedDist )
    {
        xCopyPath.exec( iCallerRunIdFrom,
                        // overlap distance from and from
                        iAllowedDist * 2, iAllowedDist * 2, //
                        iCallerRunIdTo,
                        // overlap distance between from and to
                        iAllowedDist, iAllowedDist,
                        // overlap distance to and to
                        iAllowedDist * 2, iAllowedDist * 2 );
    } // method

    inline int64_t insertCall( int64_t iSvCallerRunId, SvCall& rCall )
    {
        auto xRectangle = WKBUint64Rectangle( rCall );
        int64_t iCallId = this->insert( iSvCallerRunId, //
                                        (uint32_t)rCall.xXAxis.start( ), //
                                        (uint32_t)rCall.xYAxis.start( ), //
                                        (uint32_t)rCall.xXAxis.size( ), //
                                        (uint32_t)rCall.xYAxis.size( ), //
                                        rCall.bFromForward, //
                                        rCall.bToForward,
                                        // can deal with nullpointers
                                        makeSharedCompNucSeq( rCall.pInsertedSequence ), //
                                        rCall.pInsertedSequence == nullptr ? 0 : rCall.pInsertedSequence->length( ),
                                        (uint32_t)rCall.uiNumSuppReads, (uint32_t)rCall.uiSuppNt,
                                        (uint32_t)rCall.uiReferenceAmbiguity, rCall.iOrderID, rCall.iCtgOrderID,
                                        rCall.bMirrored, xRectangle );
        rCall.iId = iCallId;

        return iCallId;
    } // method

    inline int64_t updateCall( int64_t iSvCallerRunId, SvCall& rCall )
    {
        auto xRectangle = WKBUint64Rectangle( rCall );
        xUpdateCall.exec(
            (uint32_t)rCall.xXAxis.start( ), (uint32_t)rCall.xYAxis.start( ), (uint32_t)rCall.xXAxis.size( ),
            (uint32_t)rCall.xYAxis.size( ), rCall.bFromForward, rCall.bToForward,
            // can deal with nullpointers
            makeSharedCompNucSeq( rCall.pInsertedSequence ),
            rCall.pInsertedSequence == nullptr ? 0 : rCall.pInsertedSequence->length( ), (uint32_t)rCall.uiNumSuppReads,
            (uint32_t)rCall.uiSuppNt, (uint32_t)rCall.uiReferenceAmbiguity, rCall.iOrderID, rCall.iCtgOrderID,
            rCall.bMirrored, xRectangle, rCall.iId );
        return rCall.iId;
    } // method

    inline int64_t callArea( int64_t iCallerRunId, double dMinScore )
    {
        return xCallArea.scalar( iCallerRunId, dMinScore );
    } // method

    inline double maxScore( int64_t iCallerRunId )
    {
        return xMaxScore.scalar( iCallerRunId );
    } // method

    inline double minScore( int64_t iCallerRunId )
    {
        return xMinScore.scalar( iCallerRunId );
    } // method

    inline double maxAvSuppNt( int64_t iCallerRunId )
    {
        return xMaxAvSuppNt.scalar( iCallerRunId );
    } // method

    inline int64_t filterCallsWithHighScore( int64_t iCallerRunId, double dPercentToFilter )
    {
        if( numCalls( iCallerRunId, 0 ) == 0 )
            return 0;
        auto xTransaction = this->pConnection->sharedGuardedTrxn( );
        double dMinScore = minScore( iCallerRunId );
        double dMaxScore = maxScore( iCallerRunId );
        return xFilterCallsWithHighScore.exec( iCallerRunId,
                                               dMinScore + ( dMaxScore - dMinScore ) * ( 1 - dPercentToFilter ) );
    } // method


    using NextCallType = SQLQuery<DBCon, int64_t, bool, bool, uint32_t, uint32_t, uint32_t, uint32_t,
                                  std::shared_ptr<CompressedNucSeq>, bool, uint32_t, int64_t>;

    /** @brief returns call id, jump start pos, next context, inserted sequence, jump end position, one_sided_mate_id,
     * one sided mate do_reverse_context, last context, new contig order
     *  @details helper function for reconstructSequencedGenome */
    inline std::tuple<int64_t, uint32_t, bool, std::shared_ptr<NucSeq>, uint32_t, bool, int64_t>
    getNextCall( int64_t iOrderId, int64_t iCtgOrderId, NextCallType& xNextCall )
    {
        std::tuple<int64_t, uint32_t, bool, std::shared_ptr<NucSeq>, uint32_t, bool, int64_t> xRet;
        std::get<0>( xRet ) = -1;
        std::get<3>( xRet ) = nullptr;

        if( xNextCall.execAndFetch( iOrderId, iCtgOrderId ) )
        {
            auto xQ = xNextCall.get( );
#if 0
            std::cout << " id: " << std::get<0>( xQ )
                      << " from_forward: " << std::get<1>( xQ ) << " to_forward: " << std::get<2>( xQ )
                      << " from_pos: " << std::get<3>( xQ ) << " to_pos: " << std::get<4>( xQ )
                      << " from_size: " << std::get<5>( xQ ) << " to_size: " << std::get<6>( xQ )
                      << " inserted_sequence: "
                      << ( std::get<7>( xQ ) == nullptr ? "NULL" : std::get<7>( xQ )->pUncomNucSeq->toString( ) )
                      << " inserted_sequence_size: " << std::get<9>( xQ ) << " do_reverse: " << std::get<8>( xQ )
                      << std::endl;
#endif
            // if the call was reverted
            if( std::get<8>( xQ ) )
            {
                // reverse the call
                std::get<1>( xQ ) = !std::get<1>( xQ );
                std::get<2>( xQ ) = !std::get<2>( xQ );
                std::swap( std::get<1>( xQ ), std::get<2>( xQ ) );

                std::swap( std::get<3>( xQ ), std::get<4>( xQ ) );
                std::swap( std::get<5>( xQ ), std::get<6>( xQ ) );
            } // if

            std::get<0>( xRet ) = std::get<0>( xQ );
            // if we are in a forward context, the start of the jump is at the right of the call
            // if we are in a backward context, the start of the jump is at the left of the call
            std::get<1>( xRet ) = std::get<1>( xQ ) ? std::get<3>( xQ ) + std::get<5>( xQ ) : std::get<3>( xQ );
            // next context is simply the output of the call
            std::get<2>( xRet ) = std::get<2>( xQ );
            // check if we have an insertion
            std::get<3>( xRet ) = std::get<7>( xQ ) == nullptr ? nullptr : std::get<7>( xQ )->pUncomNucSeq;
            if( std::get<3>( xRet ) != nullptr && std::get<9>( xQ ) != std::get<3>( xRet )->length( ) )
                throw std::runtime_error( "sanity check failed: inserted_sequence is inconsistent with "
                                          "inserted_sequence_size column for call with id " +
                                          std::to_string( std::get<0>( xQ ) ) +
                                          " lengths are: " + std::to_string( std::get<9>( xQ ) ) + " and " +
                                          std::to_string( std::get<3>( xRet )->length( ) ) );

            // if we have a forward context next, the end of the jump is at the bottom of the call
            // if we have a backward context next, the end of the jump is at the top of the call
            std::get<4>( xRet ) = std::get<2>( xQ ) ? std::get<4>( xQ ) : std::get<4>( xQ ) + std::get<6>( xQ );

            // previous context
            std::get<5>( xRet ) = std::get<1>( xQ );
            // next contig order
            std::get<6>( xRet ) = std::get<10>( xQ ) + 1;

            // fetch next element to terminate query
            xNextCall.next( );
        } // if
        return xRet;
    } // method

    struct ReconstructionStats
    {
        size_t uiNumProperJumps = 0;
        size_t uiNumContradictingJumps = 0;
        size_t uiNumNonConsecutiveJumps = 0;
    };

    inline std::vector<std::tuple<std::string, std::shared_ptr<Seeds>, std::vector<std::shared_ptr<NucSeq>>>>
    callsToSeedsHelper( std::shared_ptr<Pack> pRef, bool bWithInsertions, NextCallType& xNextCall,
                        ReconstructionStats* pStats = nullptr )
    {
#if 0
        auto uiNumCalls =
            SQLQuery<DBCon, uint64_t>( this->pConnection, "SELECT COUNT(*) FROM reconstruction_table" ).scalar( );
        std::cout << "num calls: " << uiNumCalls << std::endl;
#endif

        SQLQuery<DBCon, int64_t> xCallInBetweenForw(
            this->pConnection,
            "SELECT COUNT(id) "
            "FROM sv_call_table "
            "WHERE ( "
            "   (   from_pos > ? "
            "       AND from_pos < ? "
            "       AND from_forward ) "
            "   OR "
            "   (   to_pos > ? "
            "       AND to_pos < ? "
            "       AND NOT to_forward ) "
            ") "
            "AND sv_caller_run_id IN (SELECT caller_run_id FROM reconstruction_table) "
            "LIMIT 1 " );

        SQLQuery<DBCon, int64_t> xCallInBetweenRev(
            this->pConnection,
            "SELECT COUNT(id) "
            "FROM sv_call_table "
            "WHERE ( "
            "   (   from_pos < ? "
            "       AND from_pos > ? "
            "       AND NOT from_forward ) "
            "   OR "
            "   (   to_pos < ? "
            "       AND to_pos > ? "
            "       AND to_forward ) "
            ") "
            "AND sv_caller_run_id IN (SELECT caller_run_id FROM reconstruction_table) "
            "LIMIT 1 " );
        auto fIncStatsProper = [ & ]( uint32_t uiLastPos, uint32_t uiNextPos, bool bOnFrowStrand ) {
            if( pStats != nullptr )
            {
                if( ( bOnFrowStrand ? xCallInBetweenForw : xCallInBetweenRev )
                        .scalar( uiLastPos, uiNextPos, uiLastPos, uiNextPos ) > 0 )
                    pStats->uiNumNonConsecutiveJumps++;
                else
                    pStats->uiNumProperJumps++;
            } // if
        }; // lambda

        std::vector<std::tuple<std::string, std::shared_ptr<Seeds>, std::vector<std::shared_ptr<NucSeq>>>> vRet;

        int64_t uiCtgOrder = 0;
        while( true )
        {
            auto pRet = std::make_shared<Seeds>( );
            std::vector<std::shared_ptr<NucSeq>> vInsertions;
            int64_t uiOrder = 0;
            bool uiLastWasForwardContext = true;
            uint32_t uiLastPos = 0;
            nucSeqIndex uiLastEdgeInsertionSize = 0;

            while( true )
            {
                // get the next call
                std::tuple<int64_t, uint32_t, bool, std::shared_ptr<NucSeq>, uint32_t, bool, int64_t> tNextCall;
                // search for the next call that we have not visited yet...
                metaMeasureAndLogDuration<false>(
                    "SQL", [ & ]( ) { tNextCall = this->getNextCall( uiOrder, uiCtgOrder, xNextCall ); } );
#if 0
                std::cout << "id: " << std::get<0>( tNextCall ) << " currpos: " << uiLastPos
                          << " from: " << std::get<1>( tNextCall ) << " to: "
                          << std::get<4>( tNextCall ) - pRef->startOfSequenceWithId(
                                                            pRef->uiSequenceIdForPosition( std::get<4>( tNextCall ) ) )
                          << " toContig: "
                          << pRef->nameOfSequenceWithId( pRef->uiSequenceIdForPosition( std::get<4>( tNextCall ) ) )
                          << ( std::get<2>( tNextCall ) ? " forward" : " rev-comp" ) << " inserted_seq: "
                          << ( std::get<3>( tNextCall ) == nullptr
                                   ? "nullptr"
                                   : ( std::get<3>( tNextCall )->uiSize == 0 ? "empty"
                                                                             : std::get<3>( tNextCall )->toString( ) ) )
                          << std::endl;
#endif
                // if there are no more calls or the next call starts in the next chromosome
                if( std::get<0>( tNextCall ) == -1 )
                {
                    // extract the remainder of the contig we are currently in:
                    nucSeqIndex uiSize;
                    if( uiLastWasForwardContext )
                        uiSize = pRef->endOfSequenceWithIdOrReverse( pRef->uiSequenceIdForPositionOrRev( uiLastPos ) ) -
                                 uiLastPos;
                    else
                        uiSize = uiLastPos - pRef->startOfSequenceWithIdOrReverse(
                                                 pRef->uiSequenceIdForPositionOrRev( uiLastPos ) );
                    // sanity check: contig end cannot be longer than half the contigs size
                    // just ignore ends of contigs that are too long.
                    if( uiSize <
                        pRef->lengthOfSequenceWithIdOrReverse( pRef->uiSequenceIdForPositionOrRev( uiLastPos ) ) / 2 )
                    {
                        pRet->emplace_back( pRet->empty( ) ? 0 : ( pRet->back( ).end( ) + uiLastEdgeInsertionSize ),
                                            uiSize, uiLastPos, uiLastWasForwardContext );
                        if( bWithInsertions )
                            vInsertions.push_back( std::make_shared<NucSeq>( ) ); // no inserted sequence
                        uiLastEdgeInsertionSize = 0;
                    } // if
                    // stop the loop there are no more calls.
                    break;
                } // if
                // we reach this point only if there are more calls, so tNextCall is set properly here

                // if this is the first call of the chromosome
                // then add the sequence between call start and contig start
                if( uiOrder == 0 )
                {
                    bool uiLastWasForwardContext = std::get<5>( tNextCall );
                    if( uiLastWasForwardContext )
                        uiLastPos =
                            pRef->startOfSequenceWithId( pRef->uiSequenceIdForPosition( std::get<1>( tNextCall ) ) );
                    else
                        uiLastPos =
                            pRef->endOfSequenceWithId( pRef->uiSequenceIdForPosition( std::get<1>( tNextCall ) ) );
                } // if

                // check if we can reach the start of the current call with a seed
                if( uiLastWasForwardContext && uiLastPos <= std::get<1>( tNextCall ) &&
                    !pRef->bridgingPositions( uiLastPos, std::get<1>( tNextCall ) ) )
                {
                    pRet->emplace_back( pRet->empty( ) ? uiLastEdgeInsertionSize
                                                       : ( pRet->back( ).end( ) + uiLastEdgeInsertionSize ),
                                        std::get<1>( tNextCall ) - uiLastPos + 1, uiLastPos, true );
                    fIncStatsProper( uiLastPos, std::get<1>( tNextCall ), uiLastWasForwardContext );
                } // if
                else if( !uiLastWasForwardContext && uiLastPos >= std::get<1>( tNextCall ) &&
                         !pRef->bridgingPositions( uiLastPos, std::get<1>( tNextCall ) ) )
                {
                    pRet->emplace_back( pRet->empty( ) ? uiLastEdgeInsertionSize
                                                       : ( pRet->back( ).end( ) + uiLastEdgeInsertionSize ),
                                        uiLastPos - std::get<1>( tNextCall ) + 1, uiLastPos, false );
                    fIncStatsProper( uiLastPos, std::get<1>( tNextCall ), uiLastWasForwardContext );
                } // else if
                else // cannot connect last pos and this pos
                // (add zero size seed so that we do not mess up the insertion order)
                {
                    pRet->emplace_back( pRet->empty( ) ? uiLastEdgeInsertionSize
                                                       : ( pRet->back( ).end( ) + uiLastEdgeInsertionSize ),
                                        0, 0, true );
                    // first entry (uiOrder == 0) always ends up in this else scope
                    // however this is not because it is a contradiction between traversal and graph
                    // just the starting position cannot be properly set before the first entry
                    // so don't count it as a contradiction
                    if( pStats != nullptr && uiOrder != 0 )
                        pStats->uiNumContradictingJumps++;
                } // else

                // append the skipped over sequence
                if( bWithInsertions )
                {
                    if( std::get<3>( tNextCall ) != nullptr )
                        vInsertions.push_back( std::get<3>( tNextCall ) ); // have inserted sequence
                    else
                        vInsertions.push_back( std::make_shared<NucSeq>( ) ); // no inserted sequence
                } // if

                // memorize length of the insertion (even if we do not keep the nucseq)
                // so that next seed can be placed correctly
                if( std::get<3>( tNextCall ) != nullptr )
                    uiLastEdgeInsertionSize = std::get<3>( tNextCall )->length( );
                else
                    uiLastEdgeInsertionSize = 0;

                // update orderId for next fetch as well as last context and last pos
                uiLastWasForwardContext = std::get<2>( tNextCall );
                uiLastPos = std::get<4>( tNextCall );
                uiOrder = std::get<6>( tNextCall );
            } // while

            if( uiOrder == 0 ) // check if we reached the last contig
                break;

            vRet.push_back( std::make_tuple( "chr" + std::to_string( uiCtgOrder + 1 ), pRet, vInsertions ) );
            uiCtgOrder++; // move to next contig
        } // while
        return vRet;
    }

    inline std::vector<std::tuple<std::string, std::shared_ptr<Seeds>, std::vector<std::shared_ptr<NucSeq>>>>
    callsToSeedsById( std::shared_ptr<Pack> pRef, std::vector<PriKeyDefaultType> vCallerRuns, bool bWithInsertions,
                      bool bPrintStats )
    {
        auto pTransaction = this->pConnection->uniqueGuardedTrxn( );
        auto pReconstructionTable = std::make_shared<SQLTable<DBCon, PriKeyDefaultType>>(
            this->pConnection,
            json{ { TABLE_NAME, "reconstruction_table" },
                  { CPP_EXTRA, "DROP ON DESTRUCTION" },
                  { TABLE_COLUMNS,
                    {
                        { { COLUMN_NAME, "caller_run_id" },
                          { REFERENCES, "sv_caller_run_table(id)" },
                          { CONSTRAINTS, "NOT NULL" } }
                        //
                    } } } );

        // clear table
        pReconstructionTable->deleteAllRows( );

        for( auto iCallerRun : vCallerRuns )
            pReconstructionTable->insert( iCallerRun );

        NextCallType xNextCall(
            this->pConnection,
            "SELECT id, sv_call_table.from_forward, sv_call_table.to_forward, sv_call_table.from_pos, "
            "       sv_call_table.to_pos, from_size, to_size, inserted_sequence, mirrored, "
            "       inserted_sequence_size, order_id "
            "FROM sv_call_table "
            "WHERE order_id >= ? "
            "AND ctg_order_id = ? "
            "AND sv_caller_run_id IN (SELECT caller_run_id FROM reconstruction_table) "
            "ORDER BY order_id ASC "
            "LIMIT 1 " );

        ReconstructionStats xStats;

        auto xRet = callsToSeedsHelper( pRef, bWithInsertions, xNextCall, bPrintStats ? &xStats : nullptr );
        if( bPrintStats )
        {
            std::cout << "reconstruction stats:" << std::endl;
            std::cout << "Graph is in accordance with path: " << xStats.uiNumProperJumps << std::endl;
            std::cout << "Graph contradicts path: " << xStats.uiNumContradictingJumps << std::endl;
            std::cout << "Path passes an SV entry of graph: " << xStats.uiNumNonConsecutiveJumps << std::endl;
        } // if
        return xRet;
    } // method
}; // namespace libMSV

template <typename DBCon>
using CallDescTable_t = SQLTable<DBCon,
                                 int64_t, // call_id
                                 std::string // _desc_
                                 >;
const json jCallDescTableDef = { { TABLE_NAME, "call_desc_table" },
                                 { TABLE_COLUMNS, { { { COLUMN_NAME, "call_id" } }, { { COLUMN_NAME, "_desc_" } } } } };

template <typename DBCon> class CallDescTable : public CallDescTable_t<DBCon>
{
  public:
    SQLQuery<DBCon, std::string> xGetDesc;
    CallDescTable( std::shared_ptr<DBCon> pDB )
        : CallDescTable_t<DBCon>( pDB, jCallDescTableDef ),
          xGetDesc( pDB, "SELECT _desc_ FROM call_desc_table WHERE call_id = ?" )
    {} // default constructor

    inline std::string getDesc( int64_t iId )
    {
        std::string sVal = "";
        if( xGetDesc.execAndFetch( iId ) )
            sVal = xGetDesc.getVal( );
        while( !xGetDesc.eof( ) )
            xGetDesc.next( );
        return sVal;
    } // method

    void insert_py( int64_t iId, std::string sDesc )
    {
        this->insert( iId, sDesc );
    }

    inline void genIndex( )
    {
        this->addIndex( json{ { INDEX_NAME, "call_desc_index" }, { INDEX_COLUMNS, "call_id" } } );
    }
}; // class


} // namespace libMSV
