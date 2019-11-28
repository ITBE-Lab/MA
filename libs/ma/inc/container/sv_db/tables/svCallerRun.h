/**
 * @file sequencer.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "util/sqlite3.h"

namespace libMA
{

typedef CppSQLiteExtTableWithAutomaticPrimaryKey<std::string, // name
                                                 std::string, // desc
                                                 int64_t, // timestamp
                                                 int64_t // sv_jump_run_id
                                                 >
    TP_SV_CALLER_RUN_TABLE;
class SvCallerRunTable : public TP_SV_CALLER_RUN_TABLE
{
    std::shared_ptr<CppSQLiteDBExtended> pDatabase;
    CppSQLiteExtQueryStatement<int64_t> xDelete;
    CppSQLiteExtQueryStatement<int64_t> xGetId;
    CppSQLiteExtQueryStatement<std::string, std::string, int64_t, int64_t> xGetName;
    CppSQLiteExtQueryStatement<uint32_t> xNum;
    CppSQLiteExtQueryStatement<uint32_t> xExists;
    CppSQLiteExtQueryStatement<uint32_t> xNameExists;
    CppSQLiteExtQueryStatement<int64_t> xNewestUnique;
    CppSQLiteExtStatement xInsertRow2;

  public:
    SvCallerRunTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
        : TP_SV_CALLER_RUN_TABLE(
              *pDatabase, // the database where the table resides
              "sv_caller_run_table", // name of the table in the database
              // column definitions of the table
              std::vector<std::string>{"name", "desc", "time_stamp", "sv_jump_run_id"},
              // constraints for table
              std::vector<std::string>{"FOREIGN KEY (sv_jump_run_id) REFERENCES sv_jump_run_table(id)"} ),
          pDatabase( pDatabase ),
          xDelete( *pDatabase, "DELETE FROM sv_caller_run_table WHERE name == ?" ),
          xGetId( *pDatabase, "SELECT id FROM sv_caller_run_table WHERE name == ? ORDER BY time_stamp ASC LIMIT 1" ),
          xGetName( *pDatabase,
                    "SELECT name, desc, time_stamp, sv_jump_run_id FROM sv_caller_run_table WHERE id == ?" ),
          xNum( *pDatabase, "SELECT COUNT(*) FROM sv_caller_run_table " ),
          xExists( *pDatabase, "SELECT COUNT(*) FROM sv_caller_run_table WHERE id == ?" ),
          xNameExists( *pDatabase, "SELECT COUNT(*) FROM sv_caller_run_table WHERE name == ?" ),
          xNewestUnique(
              *pDatabase,
              "SELECT id FROM sv_caller_run_table AS outer WHERE ( SELECT COUNT(*) FROM sv_caller_run_table AS "
              "inner WHERE inner.name = outer.name AND inner.time_stamp >= outer.time_stamp ) < ? "
              "AND desc = ? " ),
          xInsertRow2( *pDatabase,
                       "INSERT INTO sv_caller_run_table (id, name, desc, time_stamp, sv_jump_run_id) "
                       "VALUES (NULL, ?, ?, ?, NULL)" )
    {} // default constructor

    inline void deleteName( std::string& rS )
    {
        xDelete.bindAndExecQuery<>( rS );
        // vDump( std::cout );
    } // method

    inline int64_t getId( std::string& rS )
    {
        return xGetId.scalar( rS );
    } // method

    inline bool exists( int64_t iId )
    {
        return xExists.scalar( iId ) > 0;
    } // method

    inline bool nameExists( std::string sName )
    {
        return xNameExists.scalar( sName ) > 0;
    } // method

    inline std::string getName( int64_t iId )
    {
        return std::get<0>( xGetName.vExecuteAndReturnIterator( iId ).get( ) );
    } // method

    inline std::string getDesc( int64_t iId )
    {
        return std::get<1>( xGetName.vExecuteAndReturnIterator( iId ).get( ) );
    } // method

    inline int64_t getSvJumpRunId( int64_t iId )
    {
        return std::get<3>( xGetName.vExecuteAndReturnIterator( iId ).get( ) );
    } // method

    inline std::string getDate( int64_t iId )
    {
        auto now_c = (std::time_t)std::get<2>( xGetName.vExecuteAndReturnIterator( iId ).get( ) );
        std::stringstream ss;
#ifdef _MSC_VER
#pragma warning( suppress : 4996 ) // @todo find another way to do this
        ss << std::put_time( std::localtime( &now_c ), "%c" );
#else
        ss << std::put_time( std::localtime( &now_c ), "%c" );
#endif
        return ss.str( );
    } // method

    inline uint32_t size( )
    {
        return xNum.scalar( );
    } // method

    inline int64_t insert( std::string sName, std::string sDesc, int64_t uiJumpRunId )
    {
        if( uiJumpRunId < 0 )
        {
            this->xInsertRow2.bindAndExecute(
                sName, sDesc, (int64_t)std::chrono::system_clock::to_time_t( std::chrono::system_clock::now( ) ) );
            // get the rowid = primary key of the inserted row
            return static_cast<int64_t>( pDatabase->lastRowId( ) );
        }
        return this->xInsertRow( sName, sDesc,
                                 (int64_t)std::chrono::system_clock::to_time_t( std::chrono::system_clock::now( ) ),
                                 uiJumpRunId );
    } // method

    inline std::vector<int64_t> getNewestUnique( uint32_t uiNum, std::string sDesc )
    {
        return xNewestUnique.executeAndStoreInVector<0>( uiNum, sDesc );
    } // method
}; // class

} // namespace libMA