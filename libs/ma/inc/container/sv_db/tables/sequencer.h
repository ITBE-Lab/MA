/**
 * @file sequencer.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "common.h" // NEW DATABASE INTERFACE
#include "db_sql.h"

namespace libMA
{

typedef CppSQLiteExtTableWithAutomaticPrimaryKey<std::string // sequencer name
                                                 >
    TP_SEQUENCER_TABLE;
/**
 * @brief contains the name of a sequencer run
 * @details
 * In order to create sv jumps with different parameters for different read types,
 * each read type has to have a entry in the sequencer_table.
 * Then jumps can be computed seperated for each sequencer_table entry.
 * Eventhough the jumps are created seperated, they can all be used together in the line sweep.
 */
class SequencerTable : public TP_SEQUENCER_TABLE
{
    std::shared_ptr<CppSQLiteDBExtended> pDatabase;

  public:
    SequencerTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
        : TP_SEQUENCER_TABLE( *pDatabase, // the database where the table resides
                              "sequencer_table", // name of the table in the database
                              // column definitions of the table
                              std::vector<std::string>{ "name" },
                              // constraints for table
                              std::vector<std::string>{ "UNIQUE (name)" } ),
          pDatabase( pDatabase )
    {} // default constructor

    inline int64_t insertSequencer( std::string& sSequencerName )
    {
        return xInsertRow( sSequencerName );
    } // method
}; // class

/* NEW DATABASE INTERFACE */
/* Discuss with Markus:
 * In MySQL, the constraint UNIQUE on the column "name" leads to the error message:
 * MySQL database error: BLOB/TEXT column 'name' used in key specification without a key length
 * Explanation and possible solution(s):
 * https://stackoverflow.com/questions/1827063/mysql-error-key-specification-without-a-key-length
 */

template <typename DBCon>
using _SequencerTableType = SQLTableWithAutoPriKey<DBCon,
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
template <typename DBCon> class _SequencerTable : public _SequencerTableType<DBCon>
{
  public:
    _SequencerTable( std::shared_ptr<DBCon> pDB ) : _SequencerTableType<DBCon>( pDB, jSequencerTableDef )
    {} // default constructor

    // int64_t is the type of the primary key
    inline int64_t insertSequencer( std::string& sSequencerName )
    {
        return _SequencerTableType<DBCon>::insert( sSequencerName );
    } // method
}; // class

} // namespace libMA
