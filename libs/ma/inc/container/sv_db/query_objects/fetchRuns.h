/**
 * @file fetchRuns.h
 * @brief implements libMA::SvCallerRunsFromDb that fetches information about the sv caller runs in the DB.
 * @author Markus Schmidt
 */
#include "container/sv_db/svDb.h"

#pragma once

namespace libMA
{

/**
 * @brief fetches information about the sv caller runs in the DB.
 * @details
 * A sv caller run is mereley the creation of the calls.
 * the creation of libMA::SvJump is seperate, so that multiple caller runs can use the same jumps.
 */
class SvCallerRunsFromDb
{
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<int64_t, std::string, std::string> xQuery;
    CppSQLiteExtQueryStatement<int64_t, std::string, std::string>::Iterator xTableIterator;

  public:
    /**
     * @brief queries information about all sv caller runs.
     */
    SvCallerRunsFromDb( std::shared_ptr<SV_DB> pDb )
        : pDb( pDb ),
          xQuery( *pDb->pDatabase,
                  "SELECT id, name, desc "
                  "FROM sv_caller_run_table " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( ) )
    {} // constructor

    /// @brief return the id of the current run. undefined if eof returns true
    int64_t id( )
    {
        return std::get<0>( xTableIterator.get( ) );
    } // method

    /// @brief return the name of the current run. undefined if eof returns true
    std::string name( )
    {
        return std::get<1>( xTableIterator.get( ) );
    } // method

    /// @brief return the description of the current run. undefined if eof returns true
    std::string desc( )
    {
        return std::get<2>( xTableIterator.get( ) );
    } // method

    /// @brief advances the iterator to the next run. undefined if eof returns true
    void next( )
    {
        xTableIterator.next( );
    } // method

    /// @brief returns true if there are no more runs. undefined if eof returns true
    bool eof( )
    {
        return xTableIterator.eof( );
    } // method
}; // class

} // namespace libMA

#ifdef WITH_PYTHON
/// @brief used to expose libMA::SvCallerRunsFromDb to python
void exportRunsFromDb( py::module& rxPyModuleId );
#endif