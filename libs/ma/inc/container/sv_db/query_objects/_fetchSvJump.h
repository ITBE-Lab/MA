/**
 * @file fetchSvJump.h
 * @brief implements libMA::SortedSvJumpFromSql that fetches libMA::SvJump objects from the DB.
 * @author Markus Schmidt
 */

#pragma once
#include "db_config.h"
#ifdef USE_NEW_DB_API
#include "container/sv_db/_svDb.h" // NEW DATABASE INTERFACE

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
 * @todo this should buffer the results
 */
template <typename DBCon> class SortedSvJumpFromSql
{
    const std::shared_ptr<Presetting> pSelectedSetting;

    SQLQuery<DBCon, int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, uint32_t, int64_t, int64_t>
        xQueryStart;
    SQLQuery<DBCon, int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, uint32_t, int64_t, int64_t>
        xQueryEnd;
    // SQLQuery<DBCon, int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, uint32_t, int64_t,
    //          int64_t>::Iterator xTableIteratorStart;
    // SQLQuery<DBCon, int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, uint32_t, int64_t,
    //          int64_t>::Iterator xTableIteratorEnd;

  public:
    /// @brief fetches libMA::SvJump objects from the run with id = iSvCallerRunId sorted by their start/end positions.
    SortedSvJumpFromSql( const ParameterSetManager& rParameters, std::shared_ptr<_SV_DB<DBCon>> pDb,
                         int64_t iSvCallerRunId )
        : pSelectedSetting( rParameters.getSelected( ) ),
          xQueryStart( pDb->pDatabase,
                       "SELECT sort_pos_start, from_pos, to_pos, query_from, query_to, from_forward, to_forward, "
                       "       from_seed_start, num_supporting_nt, id, read_id "
                       "FROM sv_jump_table "
                       "WHERE sv_jump_run_id = ? "
                       "ORDER BY sort_pos_start" ),
          xQueryEnd( pDb->pDatabase,
                     "SELECT sort_pos_end, from_pos, to_pos, query_from, query_to, from_forward, to_forward, "
                     "       from_seed_start, num_supporting_nt, id, read_id "
                     "FROM sv_jump_table "
                     "WHERE sv_jump_run_id = ? "
                     "ORDER BY sort_pos_end" )
    // xTableIteratorStart( xQueryStart.vExecuteAndReturnIterator( iSvCallerRunId ) ),
    // xTableIteratorEnd( xQueryEnd.vExecuteAndReturnIterator( iSvCallerRunId ) )
    {
        xQueryStart.execAndFetch( iSvCallerRunId );
        xQueryEnd.execAndFetch( iSvCallerRunId );
    } // constructor

    /**
     * @brief fetches libMA::SvJump objects.
     * @details
     * fetches libMA::SvJump objects that:
     * - are from the run with id = iSvCallerRunId
     * - are sorted by their start/end position
     * - are within the rectangle iX,iY,uiW,uiH
     */
    SortedSvJumpFromSql( const ParameterSetManager& rParameters, std::shared_ptr<_SV_DB<DBCon>> pDb,
                         int64_t iSvCallerRunId, int64_t iX, int64_t iY, uint32_t uiW, uint32_t uiH )
        : pSelectedSetting( rParameters.getSelected( ) ),
          xQueryStart( pDb->pDatabase,
                       "SELECT sort_pos_start, from_pos, to_pos, query_from, query_to, from_forward, to_forward, "
                       "       from_seed_start, num_supporting_nt, id, read_id "
                       "FROM sv_jump_table "
                       "WHERE sv_jump_run_id = ? "
                       "AND ( (from_pos >= ? AND from_pos <= ?) OR from_pos = ? ) "
                       "AND ( (to_pos >= ? AND to_pos <= ?) OR to_pos = ? ) "
                       "ORDER BY sort_pos_start" ),
          xQueryEnd( pDb->pDatabase,
                     "SELECT sort_pos_end, from_pos, to_pos, query_from, query_to, from_forward, to_forward, "
                     "       from_seed_start, num_supporting_nt, id, read_id "
                     "FROM sv_jump_table "
                     "WHERE sv_jump_run_id = ? "
                     "AND ( (from_pos >= ? AND from_pos <= ?) OR from_pos = ? ) "
                     "AND ( (to_pos >= ? AND to_pos <= ?) OR to_pos = ? ) "
                     "ORDER BY sort_pos_end" )
    // xTableIteratorStart( xQueryStart.vExecuteAndReturnIterator(
    //     iSvCallerRunId, iX, iX + uiW, std::numeric_limits<uint32_t>::max( ), iY, iY + uiH,
    //     std::numeric_limits<uint32_t>::max( ) ) ),
    // xTableIteratorEnd( xQueryEnd.vExecuteAndReturnIterator( iSvCallerRunId, iX, iX + uiW,
    //                                                         std::numeric_limits<uint32_t>::max( ), iY, iY +
    //                                                         uiH, std::numeric_limits<uint32_t>::max( ) ) )
    {
        xQueryStart.execAndFetch( iSvCallerRunId, iX, iX + uiW, std::numeric_limits<uint32_t>::max( ), iY, iY + uiH,
                                 std::numeric_limits<uint32_t>::max( ) );
        xQueryEnd.execAndFetch( iSvCallerRunId, iX, iX + uiW, std::numeric_limits<uint32_t>::max( ), iY, iY + uiH,
                               std::numeric_limits<uint32_t>::max( ) );
    } // constructor

    /**
     * @brief fetches libMA::SvJump objects.
     * @details
     * fetches libMA::SvJump objects that:
     * - are from the run with id = iSvCallerRunId
     * - start after iS (on ref)
     * - end after iE (on ref)
     */
    SortedSvJumpFromSql( const ParameterSetManager& rParameters, std::shared_ptr<_SV_DB<DBCon>> pDb,
                         int64_t iSvCallerRunId, int64_t iS, int64_t iE )
        : pSelectedSetting( rParameters.getSelected( ) ),
          xQueryStart( pDb->pDatabase,
                       "SELECT sort_pos_start, from_pos, to_pos, query_from, query_to, from_forward, to_forward, "
                       "       from_seed_start, num_supporting_nt, id, read_id "
                       "FROM sv_jump_table "
                       "WHERE sv_jump_run_id = ? "
                       "AND sort_pos_start >= ? "
                       "AND sort_pos_start <= ? "
                       "ORDER BY sort_pos_start" ),
          xQueryEnd( pDb->pDatabase,
                     "SELECT sort_pos_end, from_pos, to_pos, query_from, query_to, from_forward, to_forward, "
                     "      from_seed_start, num_supporting_nt, id, read_id "
                     "FROM sv_jump_table "
                     "WHERE sv_jump_run_id = ? "
                     "AND sort_pos_end >= ? "
                     "AND sort_pos_end <= ? "
                     "ORDER BY sort_pos_end" )
    // xTableIteratorStart( xQueryStart.vExecuteAndReturnIterator( iSvCallerRunId, iS, iE ) ),
    // xTableIteratorEnd( xQueryEnd.vExecuteAndReturnIterator( iSvCallerRunId, iS, iE ) )
    {
        assert( iE >= iS );

        xQueryStart.execAndFetch( iSvCallerRunId, iS, iE );
        xQueryEnd.execAndFetch( iSvCallerRunId, iS, iE );
    } // constructor

    /// @brief returns whether there is another jump in the start-sorted iterator
    bool hasNextStart( )
    {
        return !xQueryStart.eof( );
    } // method

    /// @brief returns whether there is another jump in the end-sorted iterator
    bool hasNextEnd( )
    {
        return !xQueryEnd.eof( );
    } // method

    /// @brief returns whether the next start-sorted jump is smaller than the next end-sorted jump.
    bool nextStartIsSmaller( )
    {
        if( !hasNextStart( ) )
            return false;
        if( !hasNextEnd( ) )
            return true;
        auto xStartTup = xQueryStart.get( );
        auto xEndTup = xQueryEnd.get( );
        return std::get<0>( xStartTup ) <= std::get<0>( xEndTup );
    } // method

    /// @brief returns the next start-sorted jump and increases the iterator
    std::shared_ptr<SvJump> getNextStart( )
    {
        assert( hasNextStart( ) );

        auto xTup = xQueryStart.get( );
        xQueryStart.next( );
        return std::make_shared<SvJump>( pSelectedSetting, std::get<1>( xTup ), std::get<2>( xTup ),
                                         std::get<3>( xTup ), std::get<4>( xTup ), std::get<5>( xTup ),
                                         std::get<6>( xTup ), std::get<7>( xTup ), std::get<8>( xTup ),
                                         std::get<9>( xTup ), std::get<10>( xTup ) );
    } // method

    /// @brief returns the next end-sorted jump and increases the iterator
    std::shared_ptr<SvJump> getNextEnd( )
    {
        assert( hasNextEnd( ) );

        auto xTup = xQueryEnd.get( );
        xQueryEnd.next( );
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

#endif // USE_NEW_DB_API
