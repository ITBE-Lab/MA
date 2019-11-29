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

typedef CppSQLiteExtTableWithAutomaticPrimaryKey<std::string, // regex
                                                 uint32_t // state
                                                 >
    TP_SV_CALL_REG_EX_TABLE;
class SvCallRegExTable : public TP_SV_CALL_REG_EX_TABLE
{
    std::shared_ptr<CppSQLiteDBExtended> pDatabase;

  public:
    SvCallRegExTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
        : TP_SV_CALL_REG_EX_TABLE( *pDatabase, // the database where the table resides
                                   "sv_call_reg_ex_table", // name of the table in the database
                                   // column definitions of the table
                                   std::vector<std::string>{"regex", "state"} ),
          pDatabase( pDatabase )
    {} // default constructor
}; // class

} // namespace libMA