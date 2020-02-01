/**
 * @file sequencer.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "container/sv_db/tables/read.h"
#include "db_config.h"

namespace libMA
{

template <typename DBCon>
using PairedReadTableType = SQLTable<DBCon, // DB connector type
                                     int64_t, // first read (foreign key)
                                     int64_t // second read (foreign key)
                                     >;

const json jPairedReadTableDef = {{TABLE_NAME, "paired_read_table"},
                                  {TABLE_COLUMNS,
                                   {{{COLUMN_NAME, "first_read"}, {CONSTRAINTS, "REFERENCES read_table(id)"}},
                                    {{COLUMN_NAME, "second_read"}, {CONSTRAINTS, "REFERENCES read_table(id)"}}}},
                                  {PRIMARY_KEY, "first_read, second_read"}};

/**
 * @brief contains the name of a sequencer run
 * @details
 * In order to create sv jumps with different parameters for different read types,
 * each read type has to have a entry in the sequencer_table.
 * Then jumps can be computed separated for each sequencer_table entry.
 * Although the jumps are created separated, they can all be used together in the line sweep.
 */
template <typename DBCon> class PairedReadTable : public PairedReadTableType<DBCon>
{
    std::shared_ptr<ReadTable<DBCon>> pReadTable;

  public:
    /* Constructor prototype */
    PairedReadTable( std::shared_ptr<DBCon> pDB, std::shared_ptr<ReadTable<DBCon>> pReadTable )
        : PairedReadTableType<DBCon>( pDB, // the database where the table resides
                                      jPairedReadTableDef ),
          pReadTable( pReadTable )
    {} // constructor

    ~PairedReadTable( )
    {} // deconstructor

    inline std::pair<int64_t, int64_t> insertRead( int64_t uiSequencerId, std::shared_ptr<NucSeq> pReadA,
                                                   std::shared_ptr<NucSeq> pReadB )
    {
        int64_t uiReadIdA = pReadTable->insertRead( uiSequencerId, pReadA );
        int64_t uiReadIdB = pReadTable->insertRead( uiSequencerId, pReadB );
        PairedReadTable<DBCon>::insert( uiReadIdA, uiReadIdB );
        return std::make_pair( uiReadIdA, uiReadIdB );
    } // method
}; // class

} // namespace libMA
