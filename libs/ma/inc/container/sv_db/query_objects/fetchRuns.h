#include "container/sv_db/svDb.h"

#pragma once

namespace libMA
{

class SvCallerRunsFromDb
{
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<int64_t, std::string, std::string> xQuery;
    CppSQLiteExtQueryStatement<int64_t, std::string, std::string>::Iterator xTableIterator;

  public:
    SvCallerRunsFromDb( std::shared_ptr<SV_DB> pDb )
        : pDb( pDb ),
          xQuery( *pDb->pDatabase,
                  "SELECT id, name, desc "
                  "FROM sv_caller_run_table " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( ) )
    {} // constructor

    int64_t id( )
    {
        return std::get<0>( xTableIterator.get( ) );
    } // method

    std::string name( )
    {
        return std::get<1>( xTableIterator.get( ) );
    } // method

    std::string desc( )
    {
        return std::get<2>( xTableIterator.get( ) );
    } // method

    void next( )
    {
        xTableIterator.next( );
    } // method

    bool eof( )
    {
        return xTableIterator.eof( );
    } // method
}; // class

} // namespace libMA

#ifdef WITH_PYTHON
void exportRunsFromDb( py::module& rxPyModuleId );
#endif