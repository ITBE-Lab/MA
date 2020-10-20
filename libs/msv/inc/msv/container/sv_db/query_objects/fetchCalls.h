/**
 * @file fetchCalls.h
 * @brief implements libMSV::SvCallsFromDb; a class that fetches libMSV::SvCall objects from a DB.
 * @author Markus Schmidt
 */
#pragma once

#include "ms/util/parameter.h"
#include "msv/container/sv_db/tables/svCall.h"
#include "msv/container/sv_db/tables/svCallSupport.h"

namespace libMSV
{
/**
 * @brief fetches libMSV::SvCall objects from a DB.
 * @details
 * Can do several 2d range queries.
 * @todo this should become a container
 */
template <typename DBCon> class SvCallsFromDb
{
    std::shared_ptr<DBCon> pConnection;
    // the following two tables are not used. However their constructors guarantee their existence and the correctness
    // of rows
    std::shared_ptr<SvCallTable<DBCon>> pSvCallTable;
    std::shared_ptr<SvCallSupportTable<DBCon>> pSvCallSupportTable;
    // rectangle for xQuery; stays uninitialized if unused
    WKBUint64Rectangle xWkb;
    SQLQuery<typename DBCon::SlaveType, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, uint32_t,
             PriKeyDefaultType, PriKeyDefaultType>
        xQuerySupport;

    struct Helper
    {
      public:
        /// @brief only fetch calls in specific area
        bool bInArea;
        /// @brief only fetch calls that overlap/do not overlap calls from different run (ground truth run)
        bool bOverlapping;
        /// @brief true: fetch overlapping calls; false: fetch non overlapping calls
        bool bWithIntersection;
        /// @brief true: do not fetch calls that are overlapped by higher scoring calls from same run
        bool bWithSelfIntersection;
        /// @brief true: ignore ground truth run calls that are overlapped by higher scoring calls
        bool bWithOtherIntersection;
        bool bWithMinScore;
        bool bWithMaxScore;
        bool bWithMinScoreGT;
        bool bWithMaxScoreGT;
        bool bJustCount;
    }; // class

    Helper xHelper;
    std::unique_ptr<SQLQuery<DBCon, PriKeyDefaultType, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, NucSeqSql,
                             uint32_t, uint32_t>>
        pQuery = nullptr;
    std::unique_ptr<SQLQuery<DBCon, uint32_t>> pQueryCount = nullptr;

    std::string rectanglesOverlapSQL( std::string sFromTable, std::string sToTable )
    {
        // clang-format off
        // ST_DWithin returns true if the geometries are within the specified distance of one another
        // makes use of indices
        return  std::string("     AND ( (ST_DWithin(") + sToTable + ".rectangle::geometry, "
                            "                       " + sFromTable + ".rectangle::geometry, ?) "
                            "           AND " + sToTable + ".from_forward = " + sFromTable + ".from_forward "
                            "           AND " + sToTable + ".to_forward = " + sFromTable + ".to_forward) "
                            // this considers reverse complemented calls
                            // Note: we need to mirror the rectangle on the matrix diagonal
                            //       and invert the strand information
                            "       OR (ST_DWithin(" + sToTable + ".rectangle::geometry, "
                            "                      " + sFromTable + ".flipped_rectangle::geometry, ?) "
                            "           AND " + sToTable + ".from_forward != " + sFromTable + ".to_forward "
                            "           AND " + sToTable + ".to_forward != " + sFromTable + ".from_forward)"
                            "     )";
        // clang-format on
    }
    std::string selfIntersectionSQL( std::string sFromTable, std::string sToTable )
    {
        // clang-format off
        return std::string("AND NOT EXISTS( "
                        // make sure that inner_table does not overlap with any other call with higher score
                        "     SELECT " ) + sToTable + ".id "
                        "     FROM sv_call_table AS " + sToTable +
                        "     WHERE " + sToTable + ".id != " + sFromTable + ".id "
                        "     AND " + sToTable + ".score >= " + sFromTable + ".score "
                        "     AND " + sToTable + ".sv_caller_run_id = " + sFromTable + ".sv_caller_run_id "
                        // Returns true if the geometries are within the specified distance of one another
                        // makes use of indices
                        + rectanglesOverlapSQL(sFromTable, sToTable) +
                        ") ";
        // clang-format on
    }

    void initQuery( Helper xNew )
    {
        bool bInitQuery = ( xNew.bJustCount && pQueryCount == nullptr ) || //
                          ( !xNew.bJustCount && pQuery == nullptr ) || //
                          xHelper.bInArea != xNew.bInArea || //
                          xHelper.bOverlapping != xNew.bOverlapping || //
                          xHelper.bWithIntersection != xNew.bWithIntersection || //
                          xHelper.bWithSelfIntersection != xNew.bWithSelfIntersection || //
                          xHelper.bWithOtherIntersection != xNew.bWithOtherIntersection || //
                          xHelper.bWithMinScore != xNew.bWithMinScore || //
                          xHelper.bWithMaxScore != xNew.bWithMaxScore || //
                          xHelper.bWithMinScoreGT != xNew.bWithMinScoreGT || //
                          xHelper.bWithMaxScoreGT != xNew.bWithMaxScoreGT || //
                          xHelper.bJustCount != xNew.bJustCount;
        if( bInitQuery )
        {
            while( pQueryCount != nullptr && !pQueryCount->eof( ) )
                pQueryCount->next( );
            pQueryCount = nullptr;
            while( pQuery != nullptr && !pQuery->eof( ) )
                pQuery->next( );
            pQuery = nullptr;
            xHelper.bInArea = xNew.bInArea;
            xHelper.bOverlapping = xNew.bOverlapping;
            xHelper.bWithIntersection = xNew.bWithIntersection;
            xHelper.bWithSelfIntersection = xNew.bWithSelfIntersection;
            xHelper.bWithOtherIntersection = xNew.bWithOtherIntersection;
            xHelper.bWithMinScore = xNew.bWithMinScore;
            xHelper.bWithMaxScore = xNew.bWithMaxScore;
            xHelper.bWithMinScoreGT = xNew.bWithMinScoreGT;
            xHelper.bWithMaxScoreGT = xNew.bWithMaxScoreGT;
            xHelper.bJustCount = xNew.bJustCount;

            std::string sQueryText =
                std::string( "FROM sv_call_table AS inner_table "
                             "WHERE sv_caller_run_id = ? " ) +
                ( !xHelper.bInArea ? "" : "AND " ST_INTERSCTS "(rectangle, ST_GeomFromWKB(?, 0)) " ) + //
                ( !xHelper.bWithMinScore ? "" : "AND score >= ? " ) + //
                ( !xHelper.bWithMaxScore ? "" : "AND score < ? " ) + //
                ( !xHelper.bWithIntersection
                      ? ""
                      : std::string( "AND " ) + ( xHelper.bOverlapping ? "" : "NOT " ) +
                            // make sure that inner_table overlaps the outer_table:
                            "EXISTS( "
                            "     SELECT outer_table.id "
                            "     FROM sv_call_table AS outer_table "
                            "     WHERE outer_table.sv_caller_run_id = ? " +
                            ( !xHelper.bWithMinScoreGT ? "" : "AND outer_table.score >= ? " ) + //
                            ( !xHelper.bWithMaxScoreGT ? "" : "AND outer_table.score < ? " ) +
                            rectanglesOverlapSQL( "outer_table", "inner_table" ) +
                            ( !xHelper.bWithOtherIntersection ? ""
                                                              : selfIntersectionSQL( "outer_table", "outer_table2" ) ) +
                            ") " ) +
                ( !xHelper.bWithSelfIntersection ? "" : selfIntersectionSQL( "inner_table", "inner_table2" ) );

            if( xHelper.bJustCount )
                pQueryCount = std::make_unique<SQLQuery<DBCon, uint32_t>>(
                    pConnection, std::string( "SELECT COUNT(*) " ) + sQueryText );
            else
                pQuery = std::make_unique<SQLQuery<DBCon, PriKeyDefaultType, uint32_t, uint32_t, uint32_t, uint32_t,
                                                   bool, bool, NucSeqSql, uint32_t, uint32_t>>(
                    pConnection,
                    std::string( "SELECT id, from_pos, to_pos, from_size, to_size, from_forward, to_forward, "
                                 "       inserted_sequence, supporting_reads, reference_ambiguity " ) +
                        sQueryText );
        };
    } // method

    template <typename... Types> void initFetchQuery_( Helper xNew, Types... args )
    {
        initQuery( xNew );
        pQuery->execAndFetch( args... );
    } // method

  public:
    SvCallsFromDb( std::shared_ptr<DBCon> pConnection )
        : pConnection( pConnection ),
          pSvCallTable( std::make_shared<SvCallTable<DBCon>>( pConnection ) ),
          pSvCallSupportTable( std::make_shared<SvCallSupportTable<DBCon>>( pConnection ) ),
          xQuerySupport( pConnection,
                         "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, was_mirrored, "
                         "       num_supporting_nt, sv_jump_table.id, read_id "
                         "FROM sv_call_support_table "
                         "JOIN sv_jump_table ON sv_call_support_table.jump_id = sv_jump_table.id "
                         "WHERE sv_call_support_table.call_id = ? " )
    {}

    void initFetchQuery( int64_t iSvCallerIdA, int64_t iX, int64_t iY, int64_t iW, int64_t iH, int64_t iSvCallerIdB,
                         bool bOverlapping, int64_t iAllowedDist, double dMinScore, double dMaxScore )
    {
        xWkb = geom::Rectangle<nucSeqIndex>( std::max( iX, (int64_t)0 ), std::max( iY, (int64_t)0 ),
                                             std::max( iW, (int64_t)0 ), std::max( iH, (int64_t)0 ) );
        initFetchQuery_( Helper{ .bInArea = true,
                                 .bOverlapping = bOverlapping,
                                 .bWithIntersection = true,
                                 .bWithSelfIntersection = true,
                                 .bWithOtherIntersection = false,
                                 .bWithMinScore = true,
                                 .bWithMaxScore = true,
                                 .bWithMinScoreGT = false,
                                 .bWithMaxScoreGT = false,
                                 .bJustCount = false }, //
                         iSvCallerIdA, xWkb, dMinScore, dMaxScore, iSvCallerIdB, iAllowedDist, iAllowedDist,
                         iAllowedDist, iAllowedDist );
    }

    void initFetchQuery( double dMinScoreGT, double dMaxScoreGT, int64_t iSvCallerIdA, int64_t iX, int64_t iY,
                         int64_t iW, int64_t iH, int64_t iSvCallerIdB, bool bOverlapping, int64_t iAllowedDist )
    {
        xWkb = geom::Rectangle<nucSeqIndex>( std::max( iX, (int64_t)0 ), std::max( iY, (int64_t)0 ),
                                             std::max( iW, (int64_t)0 ), std::max( iH, (int64_t)0 ) );
        initFetchQuery_( Helper{ .bInArea = true,
                                 .bOverlapping = bOverlapping,
                                 .bWithIntersection = true,
                                 .bWithSelfIntersection = false,
                                 .bWithOtherIntersection = true,
                                 .bWithMinScore = false,
                                 .bWithMaxScore = false,
                                 .bWithMinScoreGT = true,
                                 .bWithMaxScoreGT = true,
                                 .bJustCount = false }, //
                         iSvCallerIdA, xWkb, iSvCallerIdB, dMinScoreGT, dMaxScoreGT, iAllowedDist, iAllowedDist,
                         iAllowedDist, iAllowedDist );
    }

    void initFetchQuery( int64_t iSvCallerIdA, int64_t iX, int64_t iY, int64_t iW, int64_t iH, int64_t iSvCallerIdB,
                         bool bOverlapping, int64_t iAllowedDist )
    {
        xWkb = geom::Rectangle<nucSeqIndex>( std::max( iX, (int64_t)0 ), std::max( iY, (int64_t)0 ),
                                             std::max( iW, (int64_t)0 ), std::max( iH, (int64_t)0 ) );
        initFetchQuery_( Helper{ .bInArea = true,
                                 .bOverlapping = bOverlapping,
                                 .bWithIntersection = true,
                                 .bWithSelfIntersection = true,
                                 .bWithOtherIntersection = false,
                                 .bWithMinScore = false,
                                 .bWithMaxScore = false,
                                 .bWithMinScoreGT = false,
                                 .bWithMaxScoreGT = false,
                                 .bJustCount = false }, //
                         iSvCallerIdA, xWkb, iSvCallerIdB, iAllowedDist, iAllowedDist, iAllowedDist, iAllowedDist );
    }

    void initFetchQuery( int64_t iSvCallerIdA, int64_t iX, int64_t iY, int64_t iW, int64_t iH, double dMinScore,
                         double dMaxScore )
    {
        xWkb = geom::Rectangle<nucSeqIndex>( std::max( iX, (int64_t)0 ), std::max( iY, (int64_t)0 ),
                                             std::max( iW, (int64_t)0 ), std::max( iH, (int64_t)0 ) );
        initFetchQuery_( Helper{ .bInArea = true,
                                 .bOverlapping = false,
                                 .bWithIntersection = false,
                                 .bWithSelfIntersection = false,
                                 .bWithOtherIntersection = false,
                                 .bWithMinScore = true,
                                 .bWithMaxScore = true,
                                 .bWithMinScoreGT = false,
                                 .bWithMaxScoreGT = false,
                                 .bJustCount = false }, //
                         iSvCallerIdA, xWkb, dMinScore, dMaxScore );
    }

    void initFetchQuery( int64_t iSvCallerIdA, int64_t iX, int64_t iY, int64_t iW, int64_t iH )
    {
        xWkb = geom::Rectangle<nucSeqIndex>( std::max( iX, (int64_t)0 ), std::max( iY, (int64_t)0 ),
                                             std::max( iW, (int64_t)0 ), std::max( iH, (int64_t)0 ) );
        initFetchQuery_( Helper{ .bInArea = true,
                                 .bOverlapping = false,
                                 .bWithIntersection = false,
                                 .bWithSelfIntersection = false,
                                 .bWithOtherIntersection = false,
                                 .bWithMinScore = false,
                                 .bWithMaxScore = false,
                                 .bWithMinScoreGT = false,
                                 .bWithMaxScoreGT = false,
                                 .bJustCount = false }, //
                         iSvCallerIdA, xWkb );
    }

    ~SvCallsFromDb( )
    {
        while( pQueryCount != nullptr && !pQueryCount->eof( ) )
            pQueryCount->next( );
        while( pQuery != nullptr && !pQuery->eof( ) )
            pQuery->next( );
    } // destructor

    /**
     * @brief fetches the next call.
     * @details
     * behaviour is undefined if libMSV::SvCallsFromDb::hasNext returns false
     */
    SvCall next( bool bWithSupport = true )
    {
        assert( pQuery != nullptr );
        auto xTup = pQuery->get( );
        SvCall xRet( std::get<1>( xTup ), // uiFromStart
                     std::get<2>( xTup ), // uiToStart
                     std::get<3>( xTup ), // uiFromSize
                     std::get<4>( xTup ), // uiToSize
                     std::get<5>( xTup ), // from_forward
                     std::get<6>( xTup ), // to_forward
                     std::get<8>( xTup ) // supporting_reads
        );
        xRet.uiReferenceAmbiguity = std::get<9>( xTup );
        xRet.pInsertedSequence = std::get<7>( xTup ).pNucSeq;
        xRet.iId = std::get<0>( xTup );

        if( bWithSupport )
        {
            xQuerySupport.execAndFetch( std::get<0>( xTup ) );
            while( !xQuerySupport.eof( ) )
            {
                auto xTup = xQuerySupport.get( );
                xRet.vSupportingJumpIds.push_back( std::get<8>( xTup ) );
                xRet.vSupportingJumps.push_back( std::make_shared<SvJump>(
                    std::get<0>( xTup ), std::get<1>( xTup ), std::get<2>( xTup ), std::get<3>( xTup ),
                    std::get<4>( xTup ), std::get<5>( xTup ), std::get<6>( xTup ), std::get<7>( xTup ),
                    std::get<8>( xTup ), std::get<9>( xTup ) ) );
                xQuerySupport.next( );
            } // while
        } // if

        pQuery->next( );
        return xRet;
    } // method

    /**
     * @brief returns true if there is another call to fetch using libMSV::SvCallsFromDb::next
     */
    bool hasNext( )
    {
        assert( pQuery != nullptr );
        return !pQuery->eof( );
    } // method

    // pair(list(tuple(x, calls with score > x, true positives with score > x)), |gt|)
    std::pair<std::vector<std::tuple<double, uint32_t, uint32_t>>, uint32_t>
    count( int64_t iSvCallerIdA, int64_t iSvCallerIdB, int64_t iAllowedDist, double dMinScore, double dMaxScore,
           double dStep )
    {
        std::vector<std::tuple<double, uint32_t, uint32_t>> vRet;
        initQuery( Helper{ false, true, true, true, false, true, true, false, false, true } );
        for( double dCurr = dMinScore; dCurr < dMaxScore; dCurr += dStep )
        {
            vRet.emplace_back( dCurr,
                               pSvCallTable->numCalls( iSvCallerIdA, dCurr ),
                               pQueryCount->scalar( iSvCallerIdA, dCurr, dMaxScore, iSvCallerIdB, iAllowedDist,
                                                    iAllowedDist, iAllowedDist, iAllowedDist ) );
        }
        return std::make_pair( vRet, pSvCallTable->numCalls( iSvCallerIdB, 0 ) );
    }
}; // namespace libMSV

} // namespace libMSV


#ifdef WITH_PYTHON
/**
 * @brief used to expose libMSV::SvCallsFromDb to python
 */
void exportCallsFromDb( libMS::SubmoduleOrganizer& xOrganizer );
#endif
