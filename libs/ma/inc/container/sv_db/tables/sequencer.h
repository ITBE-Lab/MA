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

typedef CppSQLiteExtTableWithAutomaticPrimaryKey<std::string // sequencer name
                                                 >
    TP_SEQUENCER_TABLE;
class SequencerTable : public TP_SEQUENCER_TABLE
{
    std::shared_ptr<CppSQLiteDBExtended> pDatabase;

  public:
    SequencerTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
        : TP_SEQUENCER_TABLE( *pDatabase, // the database where the table resides
                              "sequencer_table", // name of the table in the database
                              // column definitions of the table
                              std::vector<std::string>{"name"},
                              // constraints for table
                              std::vector<std::string>{"UNIQUE (name)"} ),
          pDatabase( pDatabase )
    {} // default constructor

    inline int64_t insertSequencer( std::string& sSequencerName )
    {
        return xInsertRow( sSequencerName );
    } // method
}; // class

} // namespace libMA