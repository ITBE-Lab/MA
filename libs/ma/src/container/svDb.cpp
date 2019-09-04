#include "container/svDb.h"
#include "module/combineOverlappingCalls.h"

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
        .def( "get_blur_on_overlaps_between_calls", &SV_DB::getBlurOnOverlapsBetweenCalls )
        .def( "get_num_invalid_calls", &SV_DB::getNumInvalidCalls )
        .def( "add_score_index", &SV_DB::addScoreIndex )
        .def( "get_num_calls", &SV_DB::getNumCalls )
        .def( "get_run_name", &SV_DB::getRunName )
        .def( "get_run_desc", &SV_DB::getRunDesc )
        .def( "get_run_date", &SV_DB::getRunDate )
        .def( "get_run_jump_id", &SV_DB::getRunJumpId )
        .def( "get_num_runs", &SV_DB::getNumRuns )
        .def( "get_max_score", &SV_DB::getMaxScore )
        .def( "get_min_score", &SV_DB::getMinScore )
        .def( "get_num_nts", &SV_DB::getNumNts )
        .def( "run_exists", &SV_DB::runExists )
        .def( "name_exists", &SV_DB::nameExists )
        .def( "set_num_threads", &SV_DB::setNumThreads )
        .def( "create_jump_indices", &SV_DB::createJumpIndices )
        .def( "reconstruct_sequenced_genome", &SV_DB::reconstructSequencedGenome )
        .def( "newest_unique_runs", &SV_DB::getNewestUniqueRuns )
        .def( "update_coverage", &SV_DB::updateCoverage )
        .def( "insert_sv_caller_run", &SV_DB::insertSvCallerRun )
        .def( "insert_sv_jump_run", &SV_DB::insertSvJumpRun )
        .def( "get_read", &SV_DB::getRead )
        .def( "num_jumps", &SV_DB::numJumps );

    // export the ReadInserter class
    py::class_<SV_DB::ReadInserter, std::shared_ptr<SV_DB::ReadInserter>>( rxPyModuleId, "ReadInserter" )
        .def( py::init<std::shared_ptr<SV_DB>, std::string, std::shared_ptr<Pack>>( ) )
        .def( "insert_read", &SV_DB::ReadInserter::insertRead )
        .def_readonly( "sequencer_id", &SV_DB::ReadInserter::uiSequencerId )
        .def( "insert_fasta_files", &SV_DB::ReadInserter::insertFastaFiles )
        .def( "insert_paired_fasta_files", &SV_DB::ReadInserter::insertPairedFastaFiles )
        .def( "insert_paired_read", &SV_DB::ReadInserter::insertPairedRead );

    // export the ReadContex class
    py::class_<SV_DB::SvJumpInserter::ReadContex>( rxPyModuleId, "ReadContex" )
        .def( "insert_jump", &SV_DB::SvJumpInserter::ReadContex::insertJump );

    // export the SvJumpInserter class
    py::class_<SV_DB::SvJumpInserter, std::shared_ptr<SV_DB::SvJumpInserter>>( rxPyModuleId, "SvJumpInserter" )
        .def( py::init<std::shared_ptr<SV_DB>, std::string, std::string>( ) )
        .def( "read_context", &SV_DB::SvJumpInserter::readContext )
        .def_readonly( "sv_jump_run_id", &SV_DB::SvJumpInserter::iSvJumpRunId );

    // export the ReadContex class
    // py::class_<SV_DB::SvCallInserter::CallContex>( rxPyModuleId, "CallContex" )
    //    .def( "add_support", &SV_DB::SvCallInserter::CallContex::addSupport );

    // export the SvJumpInserter class
    py::class_<SV_DB::SvCallInserter, std::shared_ptr<SV_DB::SvCallInserter>>( rxPyModuleId, "SvCallInserter" )
        .def( py::init<std::shared_ptr<SV_DB>, int64_t>( ) )
        .def( py::init<std::shared_ptr<SV_DB>, std::string, std::string, int64_t>( ) )
        .def_readonly( "sv_caller_run_id", &SV_DB::SvCallInserter::iSvCallerRunId )
        .def( "insert_call", &SV_DB::SvCallInserter::insertCall );

    // export the SortedSvJumpFromSql class
    py::class_<SortedSvJumpFromSql>( rxPyModuleId, "SortedSvJumpFromSql" )
        .def( py::init<const ParameterSetManager&, std::shared_ptr<SV_DB>, int64_t>( ) )
        .def( py::init<const ParameterSetManager&,
                       std::shared_ptr<SV_DB>,
                       int64_t,
                       uint32_t,
                       uint32_t,
                       uint32_t,
                       uint32_t>( ) )
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
        .def( py::init<const ParameterSetManager&, std::shared_ptr<SV_DB>, int64_t>( ) )
        .def( py::init<const ParameterSetManager&, std::shared_ptr<SV_DB>, int64_t, double>( ) )
        .def( py::init<const ParameterSetManager&,
                       std::shared_ptr<SV_DB>,
                       int64_t,
                       uint32_t,
                       uint32_t,
                       uint32_t,
                       uint32_t>( ) )
        .def( "next", &SvCallsFromDb::next )
        .def( "hasNext", &SvCallsFromDb::hasNext );

    // export the NucSeqFromSql classes
    exportModule<AllNucSeqFromSql, std::shared_ptr<SV_DB>, int64_t, size_t, size_t>( rxPyModuleId, "AllNucSeqFromSql" );
    exportModule<NucSeqFromSql, std::shared_ptr<SV_DB>, int64_t>( rxPyModuleId, "NucSeqFromSql" );
    exportModule<PairedNucSeqFromSql, std::shared_ptr<SV_DB>, int64_t>( rxPyModuleId, "PairedNucSeqFromSql" );
    exportModule<SvDbInserter, std::shared_ptr<SV_DB>, std::string>(
        rxPyModuleId, "SvDbInserter", []( auto&& x ) { x.def_readonly( "jump_inserter", &SvDbInserter::xInserter ); } );
    exportModule<BufferedSvDbInserter, std::shared_ptr<SV_DB>, int64_t>(
        rxPyModuleId, "BufferedSvDbInserter", []( auto&& x ) { x.def( "commit", &BufferedSvDbInserter::commit ); } );

    rxPyModuleId.def("combine_overlapping_calls", &combineOverlappingCalls);
} // function
#endif
