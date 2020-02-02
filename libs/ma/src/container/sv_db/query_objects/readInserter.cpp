#include "container/sv_db/query_objects/readInserter.h"

using namespace libMA;

#ifdef WITH_PYTHON

using DBCon = SQLDB<MySQLConDB>;

int64_t insertReads( std::vector<std::shared_ptr<NucSeq>> vReads, std::shared_ptr<SV_Schema<DBCon>> pDb,
                      std::string sName, std::shared_ptr<Pack> pRef )
{
    _ReadInserter<DBCon> xInserter( pDb, sName, pRef );

    for( auto pRead : vReads )
        xInserter.insertRead( pRead );

    return xInserter.uiSequencerId;
} // function

void exportReadInserter( py::module& rxPyModuleId )
{
    // export the ReadInserter class
    py::class_<_ReadInserter<DBCon>, std::shared_ptr<_ReadInserter<DBCon>>>( rxPyModuleId, "ReadInserter" )
        .def( py::init<std::shared_ptr<SV_Schema<DBCon>>, std::string, std::shared_ptr<Pack>>( ) )
        .def( "insert_read", &_ReadInserter<DBCon>::insertRead )
        .def_readonly( "sequencer_id", &_ReadInserter<DBCon>::uiSequencerId )
        .def( "insert_fasta_files", &_ReadInserter<DBCon>::insertFastaFiles )
        .def( "insert_paired_fasta_files", &_ReadInserter<DBCon>::insertPairedFastaFiles )
        .def( "insert_paired_read", &_ReadInserter<DBCon>::insertPairedRead );

    rxPyModuleId.def( "insert_reads", &insertReads );
} // function

#endif