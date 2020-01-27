/**
 * @file sequencer.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "db_config.h"
#ifndef USE_NEW_DB_API

#include "db_sql.h"

namespace libMA
{

typedef CppSQLiteExtTable<int64_t, // call_id (foreign key)
                          int64_t // jump_id (foreign key)
                          >
    TP_SV_CALL_SUPPORT_TABLE;
class SvCallSupportTable : public TP_SV_CALL_SUPPORT_TABLE
{
    std::shared_ptr<CppSQLiteDBExtended> pDatabase;
    CppSQLiteExtQueryStatement<int64_t> xDeleteRun;
    CppSQLiteExtStatement xDeleteCall;

  public:
    SvCallSupportTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
        : TP_SV_CALL_SUPPORT_TABLE(
              *pDatabase, // the database where the table resides
              "sv_call_support_table", // name of the table in the database
              // column definitions of the table
              std::vector<std::string>{ "call_id", "jump_id" },
              false,
              // constraints for table
              std::vector<std::string>{ "FOREIGN KEY (call_id) REFERENCES sv_call_table(id) ON DELETE CASCADE",
                                        "FOREIGN KEY (jump_id) REFERENCES sv_jump_table(id) ON DELETE CASCADE" } ),
          pDatabase( pDatabase ),
          xDeleteRun( *pDatabase, "DELETE FROM sv_call_support_table WHERE call_id IN ( SELECT id FROM "
                                  "sv_call_table WHERE sv_caller_run_id IN ( SELECT id FROM "
                                  "sv_caller_run_table WHERE name == ?))" ),
          xDeleteCall( *pDatabase,
                       "DELETE FROM sv_call_support_table "
                       "WHERE call_id = ? " )
    {
        pDatabase->execDML( "CREATE INDEX IF NOT EXISTS sv_call_support_index ON sv_call_support_table "
                            "(call_id, jump_id)" );
    } // default constructor

    inline void deleteRun( std::string& rS )
    {
        xDeleteRun.bindAndExecQuery<>( rS );
    } // method

    inline void deleteCall( int64_t iCallId )
    {
        xDeleteCall.bindAndExecute( iCallId );
    } // method

    inline void deleteCall( SvCall& rCall )
    {
        deleteCall( rCall.iId );
    } // method
}; // class

} // namespace libMA

#else

#include "common.h"

namespace libMA
{
template <typename DBCon>
using SvCallSupportTableType = SQLTable<DBCon,
                                        int64_t, // call_id (foreign key)
                                        int64_t>; // jump_id (foreign key)

const json jSvCallSupportTableDef = {
    { TABLE_NAME, "sv_call_support_table" },
    { TABLE_COLUMNS,
      {
          { { COLUMN_NAME, "call_id" } },
          { { COLUMN_NAME, "jump_id" } },
      } },
    { FOREIGN_KEY, { { COLUMN_NAME, "call_id" }, { REFERENCES, "sv_call_table(id) ON DELETE CASCADE" } } },
    { FOREIGN_KEY, { { COLUMN_NAME, "jump_id" }, { REFERENCES, "sv_jump_table(id) ON DELETE CASCADE" } } } };


template <typename DBCon> class SvCallSupportTable : public SvCallSupportTableType<DBCon>
{
    std::shared_ptr<SQLDB<DBCon>> pDatabase;
    SQLStatement<DBCon> xDeleteRun;
    SQLStatement<DBCon> xDeleteCall;

  public:
    SvCallSupportTable( std::shared_ptr<SQLDB<DBCon>> pDatabase )
        : SvCallSupportTableType<DBCon>( pDatabase, // the database where the table resides
                                         jSvCallSupportTableDef ), // table definition
          pDatabase( pDatabase ),
          xDeleteRun( pDatabase, "DELETE FROM sv_call_support_table WHERE call_id IN ( SELECT id FROM "
                                 "sv_call_table WHERE sv_caller_run_id IN ( SELECT id FROM "
                                 "sv_caller_run_table WHERE name = ?))" ),
          xDeleteCall( pDatabase,
                       "DELETE FROM sv_call_support_table "
                       "WHERE call_id = ? " )
    {
        //DEL: pDatabase->execStmt( "CREATE INDEX IF NOT EXISTS sv_call_support_index ON sv_call_support_table "
        //DEL:                      "(call_id, jump_id)" );
		this->addIndex(json{ { "INDEX_NAME", "sv_call_support_index" }, { "INDEX_COLUMNS", "call_id, jump_id" } });
    } // default constructor

    inline void deleteRun( std::string& rS )
    {
        xDeleteRun.exec( rS );
    } // method

    inline void deleteCall( int64_t iCallId )
    {
        xDeleteCall.exec( iCallId );
    } // method

    inline void deleteCall( SvCall& rCall )
    {
        deleteCall( rCall.iId );
    } // method
}; // class

} // namespace libMA

#endif
