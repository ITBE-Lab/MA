/**
 * @file sequencer.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

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
              std::vector<std::string>{"call_id", "jump_id"},
              false,
              // constraints for table
              std::vector<std::string>{"FOREIGN KEY (call_id) REFERENCES sv_call_table(id) ON DELETE CASCADE",
                                       "FOREIGN KEY (jump_id) REFERENCES sv_jump_table(id) ON DELETE CASCADE"} ),
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