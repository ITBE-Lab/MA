#include "container/sv_db/query_objects/readInserter.h"

using namespace libMA;

#ifdef WITH_PYTHON

int64_t insertReads( std::vector<std::shared_ptr<NucSeq>> vReads, std::shared_ptr<SV_DB> pDb, std::string sName,
                     std::shared_ptr<Pack> pRef )
{
    ReadInserter xInserter( pDb, sName, pRef );

    for( auto pRead : vReads )
        xInserter.insertRead( pRead );

    return xInserter.uiSequencerId;
} // function

void exportReadInserter( py::module& rxPyModuleId )
{
    // export the ReadInserter class
    py::class_<ReadInserter, std::shared_ptr<ReadInserter>>( rxPyModuleId, "ReadInserter" )
        .def( py::init<std::shared_ptr<SV_DB>, std::string, std::shared_ptr<Pack>>( ) )
        .def( "insert_read", &ReadInserter::insertRead )
        .def_readonly( "sequencer_id", &ReadInserter::uiSequencerId )
        .def( "insert_fasta_files", &ReadInserter::insertFastaFiles )
        .def( "insert_paired_fasta_files", &ReadInserter::insertPairedFastaFiles )
        .def( "insert_paired_read", &ReadInserter::insertPairedRead );

    rxPyModuleId.def( "insert_reads", &insertReads );
} // function

#endif