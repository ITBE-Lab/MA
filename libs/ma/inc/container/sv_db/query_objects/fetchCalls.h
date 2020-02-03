/**
 * @file fetchCalls.h
 * @brief implements libMA::SvCallsFromDb; a class that fetches libMA::SvCall objects from a DB.
 * @author Markus Schmidt
 */
#pragma once

#include "container/sv_db/svSchema.h"

namespace libMA
{
/**
 * @brief fetches libMA::SvCall objects from a DB.
 * @details
 * Can do several 2d range queries.
 */
template <typename DBCon> class SvCallsFromDb
{
    const std::shared_ptr<Presetting> pSelectedSetting;
    std::shared_ptr<DBCon> pConnection;
    // the following two tables are not used. However their constructors guarantee their existence and the correctness
    // of rows
    std::shared_ptr<SvCallTable<DBCon>> pSvCallTable;
    std::shared_ptr<SvCallSupportTable<DBCon>> pSvCallSupportTable;
    // rectangle for xQuery; stays uninitialized if unused
    WKBUint64Rectangle xWkb;
    SQLQuery<DBCon, int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, NucSeqSql, uint32_t, uint32_t> xQuery;
    SQLQuery<DBCon, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, uint32_t, int64_t, int64_t> xQuerySupport;

    /// @brief called from the other constructors of this class only
    SvCallsFromDb( const ParameterSetManager& rParameters, std::shared_ptr<DBCon> pConnection )
        : pSelectedSetting( rParameters.getSelected( ) ),
          pConnection( pConnection ),
          pSvCallTable( std::make_shared<SvCallTable<DBCon>>( pConnection ) ),
          pSvCallSupportTable( std::make_shared<SvCallSupportTable<DBCon>>( pConnection ) )
    {}

  public:
    /**
     * @brief fetches all calls of a specific caller id.
     */
    SvCallsFromDb( const ParameterSetManager& rParameters, std::shared_ptr<DBCon> pConnection, int64_t iSvCallerId )
        : SvCallsFromDb( rParameters, pConnection ),
          xQuery(
              pConnection,
              "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, inserted_sequence, supporting_reads, "
              "       reference_ambiguity "
              "FROM sv_call_table "
              "WHERE sv_caller_run_id = ? " ),
          xQuerySupport( pConnection,
                         "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
                         "       num_supporting_nt, sv_jump_table.id, read_id "
                         "FROM sv_call_support_table "
                         "JOIN sv_jump_table ON sv_call_support_table.jump_id = sv_jump_table.id "
                         "WHERE sv_call_support_table.call_id = ? " )
    {
        xQuery.execAndFetch( iSvCallerId );
    } // constructor

    /**
     * @brief fetches calls depending on two caller ID's
     * @details
     * this fetches all calls with run id == iSvCallerIdA that:
     *  - bOverlapping==true: overlap with at least one call from iSvCallerIdB
     *  - bOverlapping==false: overlap with no call from iSvCallerIdB
     *  - do not overlap with a call form iSvCallerIdA with higher score
     * calls that are no further than iAllowedDist from each other are considered overlapping
     */
    SvCallsFromDb( const ParameterSetManager& rParameters, std::shared_ptr<DBCon> pConnection, int64_t iSvCallerIdA,
                   int64_t iSvCallerIdB, bool bOverlapping, int64_t iAllowedDist )
        : SvCallsFromDb( rParameters, pConnection ),
          xQuery( pConnection,
                  ( std::string( "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, inserted_sequence, "
                                 "supporting_reads, "
                                 "       reference_ambiguity "
                                 "FROM sv_call_table AS inner "
                                 "WHERE sv_caller_run_id = ? "
                                 "AND " ) +
                    ( bOverlapping ? "" : "NOT " ) +
                    // make sure that inner overlaps the outer:
                    "EXISTS( "
                    "     SELECT outer.id "
                    "     FROM sv_call_table AS outer "
                    "     WHERE idx_outer.sv_caller_run_id = ? "
                    "     AND ST_Distance(outer.rectangle, inner.rectangle) <= ? "
                    "     AND outer.switch_strand = inner.switch_strand "
                    ") "
                    // make sure that inner does not overlap with any other call with higher score
                    "AND NOT EXISTS( "
                    "     SELECT inner2.id "
                    "     FROM sv_call_table AS inner2 "
                    "     WHERE idx_inner2.id != inner.id "
                    "     AND " +
                    SvCallTable<DBCon>::getSqlForCallScore( "inner2" ) +
                    " >= " + SvCallTable<DBCon>::getSqlForCallScore( "inner" ) +
                    "     AND idx_inner2.run_id_b >= inner.id " // dim 1
                    "     AND idx_inner2.run_id_a <= inner.id " // dim 1
                    "     AND ST_Distance(inner2.rectangle, inner.rectangle) <= ? "
                    "     AND inner2.switch_strand = inner.switch_strand "
                    ") " )
                      .c_str( ) ),
          xQuerySupport( pConnection,
                         "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
                         "       num_supporting_nt, sv_jump_table.id, read_id "
                         "FROM sv_call_support_table "
                         "JOIN sv_jump_table ON sv_call_support_table.jump_id = sv_jump_table.id "
                         "WHERE sv_call_support_table.call_id = ? " )
    {
        xQuery.execAndFetch( iSvCallerIdA, iSvCallerIdB, iAllowedDist );
    } // constructor

    /**
     * @brief fetches all calls of a specific caller id with a minimum score of dMinScore.
     */
    SvCallsFromDb( const ParameterSetManager& rParameters, std::shared_ptr<DBCon> pConnection, int64_t iSvCallerId,
                   double dMinScore )
        : SvCallsFromDb( rParameters, pConnection ),
          xQuery(
              pConnection,
              "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, inserted_sequence, supporting_reads, "
              "       reference_ambiguity "
              "FROM sv_call_table "
              "WHERE sv_caller_run_id = ? "
              "AND " +
                  SvCallTable<DBCon>::getSqlForCallScore( ) + " >= ? " ),
          xQuerySupport( pConnection,
                         "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
                         "       num_supporting_nt, sv_jump_table.id, read_id "
                         "FROM sv_call_support_table "
                         "JOIN sv_jump_table ON sv_call_support_table.jump_id = sv_jump_table.id "
                         "WHERE sv_call_support_table.call_id = ? " )
    {
        xQuery.execAndFetch( iSvCallerId, dMinScore );
    } // constructor

    /**
     * @brief fetches all calls of a specific caller id in the rectangle described by uiX,uiY,uiW,uiH.
     */
    SvCallsFromDb( const ParameterSetManager& rParameters, std::shared_ptr<DBCon> pConnection, int64_t iSvCallerId,
                   uint32_t uiX, uint32_t uiY, uint32_t uiW, uint32_t uiH )
        : SvCallsFromDb( rParameters, pConnection ),
          xWkb( geomUtil::Rectangle<nucSeqIndex>( uiX, uiY, uiW, uiH ) ),
          xQuery(
              pConnection,
              "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, inserted_sequence, supporting_reads, "
              "       reference_ambiguity "
              "FROM sv_call_table "
              "WHERE sv_caller_run_id = ? "
              "AND ST_Distance(rectangle, ST_PolyFromWKB(?, 0)) = 0 " ), // @todo try ST_Intersects(g1, g2) instead
          xQuerySupport( pConnection,
                         "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
                         "       num_supporting_nt, sv_jump_table.id, read_id "
                         "FROM sv_call_support_table "
                         "JOIN sv_jump_table ON sv_call_support_table.jump_id = sv_jump_table.id "
                         "WHERE sv_call_support_table.call_id = ? " )
    {
        xQuery.execAndFetch( iSvCallerId, xWkb );
    } // constructor

    /**
     * @brief fetches all calls of a specific caller id in the rectangle uiX,uiY,uiW,uiH with a score >= dMinScore.
     */
    SvCallsFromDb( const ParameterSetManager& rParameters, std::shared_ptr<DBCon> pConnection, int64_t iSvCallerId,
                   int64_t iX, int64_t iY, int64_t iW, int64_t iH, double dMinScore )
        : SvCallsFromDb( rParameters, pConnection ),
          xWkb( geomUtil::Rectangle<nucSeqIndex>( std::max( iX, (int64_t)0 ), std::max( iY, (int64_t)0 ),
                                                  std::max( iW, (int64_t)0 ), std::max( iH, (int64_t)0 ) ) ),
          xQuery( pConnection,
                  "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, inserted_sequence, "
                  "       supporting_reads, reference_ambiguity "
                  "FROM sv_call_table "
                  "WHERE sv_caller_run_id = ? "
                  "AND ST_Distance(rectangle, ST_PolyFromWKB(?, 0)) = 0 "
                  "AND " +
                      SvCallTable<DBCon>::getSqlForCallScore( ) + " >= ? " ),
          xQuerySupport( pConnection,
                         "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
                         "       num_supporting_nt, sv_jump_table.id, read_id "
                         "FROM sv_call_support_table "
                         "JOIN sv_jump_table ON sv_call_support_table.jump_id = sv_jump_table.id "
                         "WHERE sv_call_support_table.call_id = ? " )
    {
        xQuery.execAndFetch( iSvCallerId, xWkb, dMinScore );
    } // constructor

    /**
     * @brief fetches the next call.
     * @details
     * behaviour is undefined if libMA::SvCallsFromDb::hasNext returns false
     */
    SvCall next( )
    {
        auto xTup = xQuery.get( );
        SvCall xRet( std::get<1>( xTup ), // uiFromStart
                     std::get<2>( xTup ), // uiToStart
                     std::get<3>( xTup ), // uiFromSize
                     std::get<4>( xTup ), // uiToSize
                     std::get<5>( xTup ), // bSwitchStrand
                     std::get<7>( xTup ) // supporting_reads
        );
        xRet.uiReferenceAmbiguity = std::get<8>( xTup );
        xRet.pInsertedSequence = std::get<6>( xTup ).pNucSeq;
        xRet.iId = std::get<0>( xTup );

        xQuerySupport.execAndFetch( std::get<0>( xTup ) );
        while( !xQuerySupport.eof( ) )
        {
            auto xTup = xQuerySupport.get( );
            xRet.vSupportingJumpIds.push_back( std::get<8>( xTup ) );
            xRet.vSupportingJumps.push_back( std::make_shared<SvJump>(
                pSelectedSetting, std::get<0>( xTup ), std::get<1>( xTup ), std::get<2>( xTup ), std::get<3>( xTup ),
                std::get<4>( xTup ), std::get<5>( xTup ), std::get<6>( xTup ), std::get<7>( xTup ), std::get<8>( xTup ),
                std::get<9>( xTup ) ) );
            xQuerySupport.next( );
        } // while

        xQuery.next( );
        return xRet;
    } // method

    /**
     * @brief returns true if there is another call to fetch using libMA::SvCallsFromDb::next
     */
    bool hasNext( )
    {
        // return !xTableIterator.eof( );
        return !xQuery.eof( );
    } // method
}; // class

} // namespace libMA


#ifdef WITH_PYTHON
/**
 * @brief used to expose libMA::SvCallsFromDb to python
 */
void exportCallsFromDb( py::module& rxPyModuleId );
#endif
