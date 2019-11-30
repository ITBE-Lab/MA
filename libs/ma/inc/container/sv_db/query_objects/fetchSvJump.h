/**
 * @file fetchSvJump.h
 * @brief implements libMA::SortedSvJumpFromSql that fetches libMA::SvJump objects from the DB.
 * @author Markus Schmidt
 */
#include "container/sv_db/svDb.h"

#pragma once

namespace libMA
{

/**
 * @brief fetches libMA::SvJump objects from the DB.
 * @details
 * Creates two iterators:
 * - one for sv jumps sorted by their start position (on the reference)
 * - one for sv jumps sorted by their end position (on the reference)
 * The iterators can be incremented independently.
 * This is necessary for the line sweep over the SV jumps.
 */
class SortedSvJumpFromSql
{
    const std::shared_ptr<Presetting> pSelectedSetting;
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, uint32_t, int64_t,
                               int64_t>
        xQueryStart;
    CppSQLiteExtQueryStatement<int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, uint32_t, int64_t,
                               int64_t>
        xQueryEnd;
    CppSQLiteExtQueryStatement<int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, uint32_t, int64_t,
                               int64_t>::Iterator xTableIteratorStart;
    CppSQLiteExtQueryStatement<int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, uint32_t, int64_t,
                               int64_t>::Iterator xTableIteratorEnd;

  public:
    /// @brief fetches libMA::SvJump objects from the run with id = iSvCallerRunId sorted by their start/end positions.
    SortedSvJumpFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSvCallerRunId )
        : pSelectedSetting( rParameters.getSelected( ) ),
          pDb( pDb ),
          xQueryStart( *pDb->pDatabase,
                       "SELECT sort_pos_start, from_pos, to_pos, query_from, query_to, from_forward, to_forward, "
                       "       from_seed_start, num_supporting_nt, id, read_id "
                       "FROM sv_jump_table "
                       "WHERE sv_jump_run_id == ? "
                       "ORDER BY sort_pos_start" ),
          xQueryEnd( *pDb->pDatabase,
                     "SELECT sort_pos_end, from_pos, to_pos, query_from, query_to, from_forward, to_forward, "
                     "       from_seed_start, num_supporting_nt, id, read_id "
                     "FROM sv_jump_table "
                     "WHERE sv_jump_run_id == ? "
                     "ORDER BY sort_pos_end" ),
          xTableIteratorStart( xQueryStart.vExecuteAndReturnIterator( iSvCallerRunId ) ),
          xTableIteratorEnd( xQueryEnd.vExecuteAndReturnIterator( iSvCallerRunId ) )
    {} // constructor

    /**
     * @brief fetches libMA::SvJump objects.
     * @details
     * fetches libMA::SvJump objects that:
     * - are from the run with id = iSvCallerRunId
     * - are sorted by their start/end position
     * - are within the rectangle iX,iY,uiW,uiH
     */
    SortedSvJumpFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSvCallerRunId,
                         int64_t iX, int64_t iY, uint32_t uiW, uint32_t uiH )
        : pSelectedSetting( rParameters.getSelected( ) ),
          pDb( pDb ),
          xQueryStart( *pDb->pDatabase,
                       "SELECT sort_pos_start, from_pos, to_pos, query_from, query_to, from_forward, to_forward, "
                       "       from_seed_start, num_supporting_nt, id, read_id "
                       "FROM sv_jump_table "
                       "WHERE sv_jump_run_id == ? "
                       "AND ( (from_pos >= ? AND from_pos <= ?) OR from_pos == ? ) "
                       "AND ( (to_pos >= ? AND to_pos <= ?) OR to_pos == ? ) "
                       "ORDER BY sort_pos_start" ),
          xQueryEnd( *pDb->pDatabase,
                     "SELECT sort_pos_end, from_pos, to_pos, query_from, query_to, from_forward, to_forward, "
                     "       from_seed_start, num_supporting_nt, id, read_id "
                     "FROM sv_jump_table "
                     "WHERE sv_jump_run_id == ? "
                     "AND ( (from_pos >= ? AND from_pos <= ?) OR from_pos == ? ) "
                     "AND ( (to_pos >= ? AND to_pos <= ?) OR to_pos == ? ) "
                     "ORDER BY sort_pos_end" ),
          xTableIteratorStart( xQueryStart.vExecuteAndReturnIterator(
              iSvCallerRunId, iX, iX + uiW, std::numeric_limits<uint32_t>::max( ), iY, iY + uiH,
              std::numeric_limits<uint32_t>::max( ) ) ),
          xTableIteratorEnd( xQueryEnd.vExecuteAndReturnIterator( iSvCallerRunId, iX, iX + uiW,
                                                                  std::numeric_limits<uint32_t>::max( ), iY, iY + uiH,
                                                                  std::numeric_limits<uint32_t>::max( ) ) )
    {} // constructor
    /**
     * @brief fetches libMA::SvJump objects.
     * @details
     * fetches libMA::SvJump objects that:
     * - are from the run with id = iSvCallerRunId
     * - start after iS (on ref)
     * - end after iE (on ref)
     */
    SortedSvJumpFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSvCallerRunId,
                         int64_t iS, int64_t iE )
        : pSelectedSetting( rParameters.getSelected( ) ),
          pDb( pDb ),
          xQueryStart( *pDb->pDatabase,
                       "SELECT sort_pos_start, from_pos, to_pos, query_from, query_to, from_forward, to_forward, "
                       "       from_seed_start, num_supporting_nt, id, read_id "
                       "FROM sv_jump_table "
                       "WHERE sv_jump_run_id == ? "
                       "AND sort_pos_start >= ? "
                       "AND sort_pos_start <= ? "
                       "ORDER BY sort_pos_start" ),
          xQueryEnd( *pDb->pDatabase,
                     "SELECT sort_pos_end, from_pos, to_pos, query_from, query_to, from_forward, to_forward, "
                     "      from_seed_start, num_supporting_nt, id, read_id "
                     "FROM sv_jump_table "
                     "WHERE sv_jump_run_id == ? "
                     "AND sort_pos_end >= ? "
                     "AND sort_pos_end <= ? "
                     "ORDER BY sort_pos_end" ),
          xTableIteratorStart( xQueryStart.vExecuteAndReturnIterator( iSvCallerRunId, iS, iE ) ),
          xTableIteratorEnd( xQueryEnd.vExecuteAndReturnIterator( iSvCallerRunId, iS, iE ) )
    {
        assert( iE >= iS );
#if DEBUG_LEVEL > 0
#if 0
        std::cout << "SortedSvJumpFromSql::xQueryStart" << std::endl;
        xQueryStart.bindAndExplain( iSvCallerRunId, iX, iX + iW );
        std::cout << "SortedSvJumpFromSql::xQueryEnd" << std::endl;
        xQueryEnd.bindAndExplain( iSvCallerRunId, iX, iX + iW );
#endif
#endif
    } // constructor

    /// @brief returns wether there is another jump in the start-sorted iterator
    bool hasNextStart( )
    {
        return !xTableIteratorStart.eof( );
    } // method

    /// @brief returns wether there is another jump in the end-sorted iterator
    bool hasNextEnd( )
    {
        return !xTableIteratorEnd.eof( );
    } // method

    /// @brief returns wether the next start-sorted jump is smaller than the next end-sorted jump.
    bool nextStartIsSmaller( )
    {
        if( !hasNextStart( ) )
            return false;
        if( !hasNextEnd( ) )
            return true;
        auto xStartTup = xTableIteratorStart.get( );
        auto xEndTup = xTableIteratorEnd.get( );
        return std::get<0>( xStartTup ) <= std::get<0>( xEndTup );
    } // method

    /// @brief returns the next start-sorted jump and increases the iterator
    std::shared_ptr<SvJump> getNextStart( )
    {
        assert( hasNextStart( ) );

        auto xTup = xTableIteratorStart.get( );
        xTableIteratorStart.next( );
        return std::make_shared<SvJump>( pSelectedSetting, std::get<1>( xTup ), std::get<2>( xTup ),
                                         std::get<3>( xTup ), std::get<4>( xTup ), std::get<5>( xTup ),
                                         std::get<6>( xTup ), std::get<7>( xTup ), std::get<8>( xTup ),
                                         std::get<9>( xTup ), std::get<10>( xTup ) );
    } // method

    /// @brief returns the next end-sorted jump and increases the iterator
    std::shared_ptr<SvJump> getNextEnd( )
    {
        assert( hasNextEnd( ) );

        auto xTup = xTableIteratorEnd.get( );
        xTableIteratorEnd.next( );
        return std::make_shared<SvJump>( pSelectedSetting, std::get<1>( xTup ), std::get<2>( xTup ),
                                         std::get<3>( xTup ), std::get<4>( xTup ), std::get<5>( xTup ),
                                         std::get<6>( xTup ), std::get<7>( xTup ), std::get<8>( xTup ),
                                         std::get<9>( xTup ), std::get<10>( xTup ) );
    } // method

}; // class

} // namespace libMA

#ifdef WITH_PYTHON
/// @brief used to expose libMA::SortedSvJumpFromSql to python
void exportSvJump( py::module& rxPyModuleId );
#endif