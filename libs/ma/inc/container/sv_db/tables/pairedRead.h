/**
 * @file sequencer.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "container/sv_db/tables/read.h"
#include "db_sql.h"

namespace libMA
{

typedef CppSQLiteExtTable<int64_t, // first read (foreign key)
                          int64_t // second read (foreign key)
                          >
    TP_PAIRED_READ_TABLE;
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