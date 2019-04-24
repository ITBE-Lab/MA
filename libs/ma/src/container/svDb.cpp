#include "container/svDb.h"

using namespace libMA;

SV_DB::PairedReadTable::PairedReadTable(
    std::shared_ptr<CppSQLiteDBExtended> pDatabase, std::shared_ptr<ReadTable> pReadTable )
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

#ifdef WITH_PYTHON

void exportSoCDbWriter( py::module& rxPyModuleId )
{

    // export the SV_DB class
    py::class_<SV_DB, std::shared_ptr<SV_DB>>( rxPyModuleId, "SV_DB" ) //
        .def( py::init<std::string, std::string>( ) )
        .def("clear_calls_table" , &SV_DB::clearCallsTable);

    // export the ReadInserter class
    py::class_<SV_DB::ReadInserter, std::shared_ptr<SV_DB::ReadInserter>>( rxPyModuleId, "ReadInserter" )
        .def( py::init<std::shared_ptr<SV_DB>, std::string>( ) )
        .def( "insert_read", &SV_DB::ReadInserter::insertRead )
        .def( "insert_paired_read", &SV_DB::ReadInserter::insertPairedRead );

    // export the ReadContex class
    py::class_<SV_DB::SvJumpInserter::ReadContex>( rxPyModuleId, "ReadContex" )
        .def( "insert_jump", &SV_DB::SvJumpInserter::ReadContex::insertJump );

    // export the SvJumpInserter class
    py::class_<SV_DB::SvJumpInserter, std::shared_ptr<SV_DB::SvJumpInserter>>( rxPyModuleId, "SvJumpInserter" )
        .def( py::init<std::shared_ptr<SV_DB>, std::string, std::string>( ) )
        .def( "insert_read", &SV_DB::SvJumpInserter::insertRead )
        .def( "read_context", &SV_DB::SvJumpInserter::readContext )
        .def_readonly( "sv_caller_run_id", &SV_DB::SvJumpInserter::iSvCallerRunId );

    // export the NucSeqFromSql classes
    exportModule<AllNucSeqFromSql, std::shared_ptr<SV_DB>>( rxPyModuleId, "AllNucSeqFromSql" );
    exportModule<NucSeqFromSql, std::shared_ptr<SV_DB>>( rxPyModuleId, "NucSeqFromSql" );
    exportModule<PairedNucSeqFromSql, std::shared_ptr<SV_DB>>( rxPyModuleId, "PairedNucSeqFromSql" );
} // function
#endif
