/**
 * @file fetchRuns.h
 * @brief implements libMA::SvCallerRunsFromDb that fetches information about the sv caller runs in the DB.
 * @author Markus Schmidt
 */
#pragma once

#include "container/sv_db/svSchema.h"
#include "db_config.h"

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
    std::shared_ptr<SV_Schema<DBCon>> pDb;
    SQLQuery<DBCon, int64_t, std::string, std::string> xQuery;
    // DELETED: SQLQuery<DBCon, int64_t, std::string, std::string>::Iterator xTableIterator;

  public:
    /**
     * @brief queries information about all sv caller runs.
     */
    SvCallerRunsFromDb( std::shared_ptr<SV_Schema<DBCon>> pDb )
        : pDb( pDb ),
          xQuery( pDb->pDatabase,
                  "SELECT id, name, _desc_ "
                  "FROM sv_caller_run_table" )
    // DELETED: xTableIterator(xQuery.vExecuteAndReturnIterator())
    {
        xQuery.execAndFetch( );
        int64_t iId = std::get<0>( xQuery.get( ) );
        std::cout << "\nValue for iId after execution is: " << iId << std::endl;
    } // constructor

    /// @brief return the id of the current run. undefined if eof returns true
    int64_t id( )
    {
        int64_t iId = std::get<0>( xQuery.get( ) );
        std::cout << "\nValue for iId is: " << iId << std::endl;
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