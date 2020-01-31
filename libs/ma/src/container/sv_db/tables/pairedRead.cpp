#include "db_config.h"
#ifndef USE_NEW_DB_API
#include "container/sv_db/svDb.h"


#include "container/sv_db/tables/pairedRead.h"

using namespace libMA;

PairedReadTable::PairedReadTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase,
                                  std::shared_ptr<ReadTable>
                                      pReadTable )
    : TP_PAIRED_READ_TABLE( *pDatabase, // the database where the table resides
                            "paired_read_table", // name of the table in the database
                            // column definitions of the table
                            std::vector<std::string>{"first_read", "second_read"},
                            false, // do not generate automatic primary key...
                            // constraints for table
                            std::vector<std::string>{"FOREIGN KEY (first_read) REFERENCES read_table(id)", //
                                                     "FOREIGN KEY (second_read) REFERENCES read_table(id)", //
                                                     "PRIMARY KEY (first_read, second_read)"} ),
      pDatabase( pDatabase ),
      pReadTable( pReadTable )
{} // default constructor

#else
#include "container/sv_db/_svDb.h" // NEW DATABASE INTERFACE
#endif