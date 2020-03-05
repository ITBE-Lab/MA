/**
 * @file sequencer.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "sql_api.h" // NEW DATABASE INTERFACE

namespace libMSV
{

/* NEW DATABASE INTERFACE */
/* Discuss with Markus:
 * In MySQL, the constraint UNIQUE on the column "name" leads to the error message:
 * MySQL database error: BLOB/TEXT column 'name' used in key specification without a key length
 * Explanation and possible solution(s):
 * https://stackoverflow.com/questions/1827063/mysql-error-key-specification-without-a-key-length
 */

template <typename DBCon>
using SequencerTableType = SQLTableWithAutoPriKey<DBCon,
                                                   std::string // sequencer name
                                                   >;
const json jSequencerTableDef = { { TABLE_NAME, "sequencer_table" },
                                  { TABLE_COLUMNS, { { { COLUMN_NAME, "name" } /*, {CONSTRAINTS, "UNIQUE"} */ } } } };

/**
 * @brief contains the name of a sequencer run
 * @details
 * In order to create sv jumps with different parameters for different read types,
 * each read type has to have a entry in the sequencer_table.
 * Then jumps can be computed separated for each sequencer_table entry.
 * Although the jumps are created separated, they can all be used together in the line sweep.
 */
template <typename DBCon> class SequencerTable : public SequencerTableType<DBCon>
{
  public:
    SequencerTable( std::shared_ptr<DBCon> pDB ) : SequencerTableType<DBCon>( pDB, jSequencerTableDef )
    {} // default constructor
}; // class

} // namespace libMSV
