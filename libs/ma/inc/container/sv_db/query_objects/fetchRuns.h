/**
 * @file fetchRuns.h
 * @brief implements libMA::SvCallerRunsFromDb that fetches information about the sv caller runs in the DB.
 * @author Markus Schmidt
 */
#pragma once

#include "container/sv_db/svSchema.h"

namespace libMA
{

/**
 * @brief fetches information about the sv caller runs in the DB.
 * @details
 * A sv caller run is merely the creation of the calls.
 * the creation of libMA::SvJump is separate, so that multiple caller runs can use the same jumps.
 */
template <typename DBCon> class SvCallerRunsFromDb
{
    std::shared_ptr<DBCon> pConnection;
    // table object is not used. However its constructor guarantees its existence and the correctness of rows
    std::shared_ptr<SvCallerRunTable<DBCon>> pSvCallerRunTable;
    SQLQuery<DBCon, int64_t, std::string, std::string> xQuery;

  public:
    /**
     * @brief queries information about all sv caller runs.
     */
    SvCallerRunsFromDb( std::shared_ptr<DBCon> pConnection )
        : pConnection( pConnection ),
          pSvCallerRunTable( std::make_shared<SvCallerRunTable<DBCon>>( pConnection ) ),
          xQuery( pConnection,
                  "SELECT id, name, _desc_ "
                  "FROM sv_caller_run_table" )
    {
        xQuery.execAndFetch( );
    } // constructor

    /// @brief return the id of the current run. undefined if eof returns true
    int64_t id( )
    {
        int64_t iId = std::get<0>( xQuery.get( ) );
        return iId;
    } // method

    /// @brief return the name of the current run. undefined if eof returns true
    std::string name( )
    {
        return std::get<1>( xQuery.get( ) );
    } // method

    /// @brief return the description of the current run. undefined if eof returns true
    std::string desc( )
    {
        return std::get<2>( xQuery.get( ) );
    } // method

    /// @brief advances the iterator to the next run. undefined if eof returns true
    void next( )
    {
        xQuery.next( );
    } // method

    /// @brief returns true if there are no more runs. undefined if eof returns true
    bool eof( )
    {
        return xQuery.eof( );
    } // method
}; // class

} // namespace libMA

#ifdef WITH_PYTHON
/// @brief used to expose libMA::SvCallerRunsFromDb to python
void exportRunsFromDb( py::module& rxPyModuleId );
#endif