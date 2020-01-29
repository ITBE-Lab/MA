/**
 * @file sequencer.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "container/sv_db/tables/read.h"

#include "db_config.h"
#ifndef USE_NEW_DB_API

#include "db_sql.h"
namespace libMA
{

typedef CppSQLiteExtTable<int64_t, // first read (foreign key)
                          int64_t // second read (foreign key)
                          >
    TP_PAIRED_READ_TABLE;
/**
 * @brief this table connects reads by their ids in order to indicate paired reads
 */
class PairedReadTable : public TP_PAIRED_READ_TABLE
{
    std::shared_ptr<CppSQLiteDBExtended> pDatabase;
    std::shared_ptr<ReadTable> pReadTable;

  public:
    PairedReadTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase, std::shared_ptr<ReadTable> pReadTable );

    ~PairedReadTable( )
    {} // deconstructor

    inline std::pair<int64_t, int64_t> insertRead( int64_t uiSequencerId, std::shared_ptr<NucSeq> pReadA,
                                                   std::shared_ptr<NucSeq> pReadB )
    {
        int64_t uiReadIdA = pReadTable->insertRead( uiSequencerId, pReadA );
        int64_t uiReadIdB = pReadTable->insertRead( uiSequencerId, pReadB );
        xInsertRow( uiReadIdA, uiReadIdB );
        return std::make_pair( uiReadIdA, uiReadIdB );
    } // method
}; // class

} // namespace libMA
#else

namespace libMA
{
/* NEW DATABASE INTERFACE */

template <typename DBCon>
using PairedReadTableType = SQLTable<DBCon, // DB connector type
                                     int64_t, // first read (foreign key)
                                     int64_t // second read (foreign key)
                                     >;

const json jPairedReadTableDef = {
    { TABLE_NAME, "paired_read_table" },
    { TABLE_COLUMNS,
      { { { COLUMN_NAME, "first_read" }, { CONSTRAINTS, "REFERENCES read_table(id)" } },
        { { COLUMN_NAME, "second_read" }, { CONSTRAINTS, "REFERENCES read_table(id)" } } } },
    { PRIMARY_KEY, "first_read, second_read" } };

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
    std::shared_ptr<_ReadTable<DBCon>> pReadTable;

  public:
    /* Constructor prototype */
    PairedReadTable( std::shared_ptr<SQLDB<DBCon>> pDB, std::shared_ptr<_ReadTable<DBCon>> pReadTable )
        : PairedReadTableType<DBCon>( pDB, // the database where the table resides
                                      jPairedReadTableDef )
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
#endif
