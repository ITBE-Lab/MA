#include "container/sv_db/svDb.h"

#pragma once

namespace libMA
{

class SvCallsFromDb
{
    const std::shared_ptr<Presetting> pSelectedSetting;
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, NucSeqSql, uint32_t, uint32_t>
        xQuery;
    CppSQLiteExtQueryStatement<uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, uint32_t, int64_t, int64_t>
        xQuerySupport;
    CppSQLiteExtQueryStatement<int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, NucSeqSql, uint32_t,
                               uint32_t>::Iterator xTableIterator;

  public:
    SvCallsFromDb( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSvCallerId )
        : pSelectedSetting( rParameters.getSelected( ) ),
          pDb( pDb ),
          xQuery( *pDb->pDatabase,
                  "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, inserted_sequence, supporting_nt, "
                  "       coverage "
                  "FROM sv_call_table "
                  "WHERE sv_caller_run_id == ? " ),
          xQuerySupport( *pDb->pDatabase,
                         "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
                         "num_supporting_nt, sv_jump_table.id, read_id "
                         "FROM sv_call_support_table "
                         "JOIN sv_jump_table ON sv_call_support_table.jump_id == sv_jump_table.id "
                         "WHERE sv_call_support_table.call_id == ? " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( iSvCallerId ) )
    {} // constructor

    // fetch overlapping or non overlapping calls:
    SvCallsFromDb( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSvCallerIdA,
                   int64_t iSvCallerIdB, bool bOverlapping, int64_t iAllowedDist )
        : pSelectedSetting( rParameters.getSelected( ) ),
          pDb( pDb ),
          xQuery(
              *pDb->pDatabase,
              ( std::string(
                    "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, inserted_sequence, supporting_nt, "
                    "       coverage "
                    "FROM sv_call_table AS inner "
                    "WHERE sv_caller_run_id == ? "
                    "AND " ) +
                ( bOverlapping ? "" : "NOT " ) +
                // make sure that inner overlaps the outer:
                "EXISTS( "
                "     SELECT outer.id "
                "     FROM sv_call_table AS outer, sv_call_r_tree AS idx_outer "
                "     WHERE outer.id == idx_outer.id "
                "     AND idx_outer.run_id_b >= ? " // dim 1
                "     AND idx_outer.run_id_a <= ? " // dim 1
                "     AND idx_outer.maxX >= inner.from_pos - ? " // dim 2
                "     AND idx_outer.minX <= inner.from_pos + inner.from_size + ? " // dim 2
                "     AND idx_outer.maxY >= inner.to_pos - ? " // dim 3
                "     AND idx_outer.minY <= inner.to_pos + inner.to_size + ? " // dim 3
                "     AND outer.switch_strand == inner.switch_strand "
                ") "
                // make sure that inner does not overlap with any other call with higher score
                "AND NOT EXISTS( "
                "     SELECT inner2.id "
                "     FROM sv_call_table AS inner2, sv_call_r_tree AS idx_inner2 "
                "     WHERE inner2.id == idx_inner2.id "
                "     AND idx_inner2.id != inner.id "
                "     AND (inner2.supporting_nt*1.0)/inner2.coverage >= (inner.supporting_nt*1.0)/inner.coverage "
                "     AND idx_inner2.run_id_b >= inner.id " // dim 1
                "     AND idx_inner2.run_id_a <= inner.id " // dim 1
                "     AND idx_inner2.maxX >= inner.from_pos - ? " // dim 2
                "     AND idx_inner2.minX <= inner.from_pos + inner.from_size + ? " // dim 2
                "     AND idx_inner2.maxY >= inner.to_pos - ? " // dim 3
                "     AND idx_inner2.minY <= inner.to_pos + inner.to_size + ? " // dim 3
                "     AND inner2.switch_strand == inner.switch_strand "
                ") " )
                  .c_str( ) ),
          xQuerySupport( *pDb->pDatabase,
                         "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
                         "num_supporting_nt, sv_jump_table.id, read_id "
                         "FROM sv_call_support_table "
                         "JOIN sv_jump_table ON sv_call_support_table.jump_id == sv_jump_table.id "
                         "WHERE sv_call_support_table.call_id == ? " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( iSvCallerIdA, iSvCallerIdB, iSvCallerIdB, iAllowedDist,
                                                            iAllowedDist, iAllowedDist, iAllowedDist, iAllowedDist,
                                                            iAllowedDist, iAllowedDist, iAllowedDist ) )
    {} // constructor

    SvCallsFromDb( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSvCallerId,
                   double dMinScore )
        : pSelectedSetting( rParameters.getSelected( ) ),
          pDb( pDb ),
          xQuery( *pDb->pDatabase,
                  "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, inserted_sequence, supporting_nt, "
                  "       coverage "
                  "FROM sv_call_table "
                  "WHERE sv_caller_run_id == ? "
                  "AND (supporting_nt*1.0)/coverage >= ? " ),
          xQuerySupport( *pDb->pDatabase,
                         "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
                         "num_supporting_nt, sv_jump_table.id, read_id "
                         "FROM sv_call_support_table "
                         "JOIN sv_jump_table ON sv_call_support_table.jump_id == sv_jump_table.id "
                         "WHERE sv_call_support_table.call_id == ? " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( iSvCallerId, dMinScore ) )
    {} // constructor

    SvCallsFromDb( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSvCallerId,
                   uint32_t uiX, uint32_t uiY, uint32_t uiW, uint32_t uiH )
        : pSelectedSetting( rParameters.getSelected( ) ),
          pDb( pDb ),
          xQuery( *pDb->pDatabase,
                  "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, inserted_sequence, supporting_nt, "
                  "       coverage "
                  "FROM sv_call_table "
                  "WHERE sv_caller_run_id == ? "
                  "AND from_pos + from_size >= ? "
                  "AND to_pos + to_size >= ? "
                  "AND from_pos <= ? "
                  "AND to_pos <= ? " ),
          xQuerySupport( *pDb->pDatabase,
                         "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
                         "num_supporting_nt, sv_jump_table.id, read_id "
                         "FROM sv_call_support_table "
                         "JOIN sv_jump_table ON sv_call_support_table.jump_id == sv_jump_table.id "
                         "WHERE sv_call_support_table.call_id == ? " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( iSvCallerId, uiX, uiY, uiX + uiW, uiY + uiH ) )
    {} // constructor

    SvCallsFromDb( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSvCallerId, int64_t iX,
                   int64_t iY, int64_t iW, int64_t iH, double dMinScore )
        : pSelectedSetting( rParameters.getSelected( ) ),
          pDb( pDb ),
          xQuery( *pDb->pDatabase,
                  "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, inserted_sequence, supporting_nt, "
                  "       coverage "
                  "FROM sv_call_table "
                  "WHERE sv_caller_run_id == ? "
                  "AND from_pos + from_size >= ? "
                  "AND to_pos + to_size >= ? "
                  "AND from_pos <= ? "
                  "AND to_pos <= ? "
                  "AND (supporting_nt*1.0)/coverage >= ? " ),
          xQuerySupport( *pDb->pDatabase,
                         "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
                         "num_supporting_nt, sv_jump_table.id, read_id "
                         "FROM sv_call_support_table "
                         "JOIN sv_jump_table ON sv_call_support_table.jump_id == sv_jump_table.id "
                         "WHERE sv_call_support_table.call_id == ? " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( iSvCallerId, iX, iY, iX + iW, iY + iH, dMinScore ) )
    {} // constructor

    SvCall next( )
    {
        auto xTup = xTableIterator.get( );
        SvCall xRet( std::get<1>( xTup ), // uiFromStart
                     std::get<2>( xTup ), // uiToStart
                     std::get<3>( xTup ), // uiFromSize
                     std::get<4>( xTup ), // uiToSize
                     std::get<5>( xTup ), // bSwitchStrand
                     std::get<7>( xTup ) // num_supporting_nt
        );
        xRet.uiCoverage = std::get<8>( xTup );
        xRet.pInsertedSequence = std::get<6>( xTup ).pNucSeq;
        xRet.iId = std::get<0>( xTup );
        auto xSupportIterator( xQuerySupport.vExecuteAndReturnIterator( std::get<0>( xTup ) ) );
        while( !xSupportIterator.eof( ) )
        {
            auto xTup = xSupportIterator.get( );
            xRet.vSupportingJumpIds.push_back( std::get<7>( xTup ) );
            xRet.vSupportingJumps.push_back( std::make_shared<SvJump>(
                pSelectedSetting, std::get<0>( xTup ), std::get<1>( xTup ), std::get<2>( xTup ), std::get<3>( xTup ),
                std::get<4>( xTup ), std::get<5>( xTup ), std::get<6>( xTup ), std::get<7>( xTup ), std::get<8>( xTup ),
                std::get<9>( xTup ) ) );
            xSupportIterator.next( );
        } // while
        xTableIterator.next( );
        return xRet;
    } // method

    bool hasNext( )
    {
        return !xTableIterator.eof( );
    } // method
}; // class

} // namespace libMA

#ifdef WITH_PYTHON
void exportCallsFromDb( py::module& rxPyModuleId );
#endif