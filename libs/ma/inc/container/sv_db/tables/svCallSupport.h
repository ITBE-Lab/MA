/**
 * @file sequencer.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "common.h"
#include "container/svJump.h"

namespace libMA
{
template <typename DBCon>
using SvCallSupportTableType = SQLTable<DBCon,
                                        PriKeyDefaultType, // call_id (foreign key)
                                        PriKeyDefaultType>; // jump_id (foreign key)

const json jSvCallSupportTableDef = {
    {TABLE_NAME, "sv_call_support_table"},
    {TABLE_COLUMNS,
     {
         {{COLUMN_NAME, "call_id"}},
         {{COLUMN_NAME, "jump_id"}},
     }},
    {FOREIGN_KEY, {{COLUMN_NAME, "call_id"}, {REFERENCES, "sv_call_table(id) ON DELETE CASCADE"}}},
    {FOREIGN_KEY, {{COLUMN_NAME, "jump_id"}, {REFERENCES, "sv_jump_table(id) ON DELETE CASCADE"}}}};


template <typename DBCon> class SvCallSupportTable : public SvCallSupportTableType<DBCon>
{
    std::shared_ptr<DBCon> pDatabase;
    SQLStatement<DBCon> xDeleteRun;
    SQLStatement<DBCon> xDeleteCall;

  public:
#if DEBUG_LEVEL == 0 // @todo discuss with arne
    // increase the bulk inserter size on this table
    using uiBulkInsertSize = std::integral_constant<size_t, 5000>;
#endif

    SvCallSupportTable( std::shared_ptr<DBCon> pDatabase )
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
        // @todo reenable via some function
        // this->addIndex( json{{"INDEX_NAME", "sv_call_support_index"}, {"INDEX_COLUMNS", "call_id, jump_id"}} );
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
