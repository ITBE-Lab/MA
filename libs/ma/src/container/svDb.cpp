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
        .def( "get_run_id", &SV_DB::getRunId )
        .def( "get_call_area", &SV_DB::getCallArea )
        .def( "get_num_overlaps_between_calls", &SV_DB::getNumOverlapsBetweenCalls )
        .def( "get_num_calls", &SV_DB::getNumCalls )
        .def( "get_run_name", &SV_DB::getRunName )
        .def( "get_run_desc", &SV_DB::getRunDesc )
        .def( "get_run_date", &SV_DB::getRunDate )
        .def( "get_num_runs", &SV_DB::getNumRuns )
        .def( "run_exists", &SV_DB::runExists )
        .def( "clear_calls_table", &SV_DB::clearCallsTable )
        .def( "clear_calls_table_for_caller", &SV_DB::clearCallsTableForCaller )
        .def( "set_num_threads", &SV_DB::setNumThreads )
        .def( "create_caller_indices", &SV_DB::createCallerIndices )
        .def( "drop_caller_indices", &SV_DB::dropCallerIndices )
        .def( "create_sequencer_indices", &SV_DB::createSequencerIndices )
        .def( "num_jumps", &SV_DB::numJumps );

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

    // export the ReadContex class
    //py::class_<SV_DB::SvCallInserter::CallContex>( rxPyModuleId, "CallContex" )
    //    .def( "add_support", &SV_DB::SvCallInserter::CallContex::addSupport );

    // export the SvJumpInserter class
    py::class_<SV_DB::SvCallInserter, std::shared_ptr<SV_DB::SvCallInserter>>( rxPyModuleId, "SvCallInserter" )
        .def( py::init<std::shared_ptr<SV_DB>, int64_t>( ) )
        .def( "insert_call", &SV_DB::SvCallInserter::insertCall );

    // export the SortedSvJumpFromSql class
    py::class_<SortedSvJumpFromSql>( rxPyModuleId, "SortedSvJumpFromSql" )
        .def( py::init<std::shared_ptr<SV_DB>, int64_t>( ) )
        .def( "has_next_start", &SortedSvJumpFromSql::hasNextStart )
        .def( "has_next_end", &SortedSvJumpFromSql::hasNextEnd )
        .def( "next_start_is_smaller", &SortedSvJumpFromSql::nextStartIsSmaller )
        .def( "get_next_start", &SortedSvJumpFromSql::getNextStart )
        .def( "get_next_end", &SortedSvJumpFromSql::getNextEnd );

    // export the SvCallerRunsFromDb class
    py::class_<SvCallerRunsFromDb>( rxPyModuleId, "SvCallerRunsFromDb" )
        .def( py::init<std::shared_ptr<SV_DB>>( ) )
        .def( "id", &SvCallerRunsFromDb::id )
        .def( "name", &SvCallerRunsFromDb::name )
        .def( "desc", &SvCallerRunsFromDb::desc )
        .def( "eof", &SvCallerRunsFromDb::eof )
        .def( "next", &SvCallerRunsFromDb::next );

    // export the SvCallsFromDb class
    py::class_<SvCallsFromDb>( rxPyModuleId, "SvCallsFromDb" )
        .def( py::init<std::shared_ptr<SV_DB>, int64_t>( ) )
        .def( "next", &SvCallsFromDb::next )
        .def( "hasNext", &SvCallsFromDb::hasNext );

    // export the NucSeqFromSql classes
    exportModule<AllNucSeqFromSql, std::shared_ptr<SV_DB>>( rxPyModuleId, "AllNucSeqFromSql" );
    exportModule<NucSeqFromSql, std::shared_ptr<SV_DB>>( rxPyModuleId, "NucSeqFromSql" );
    exportModule<PairedNucSeqFromSql, std::shared_ptr<SV_DB>>( rxPyModuleId, "PairedNucSeqFromSql" );
    exportModule<SvDbInserter, std::shared_ptr<SV_DB>, std::string>(
        rxPyModuleId, "SvDbInserter", []( auto&& x ) { x.def_readonly( "jump_inserter", &SvDbInserter::xInserter ); } );
} // function
#endif
