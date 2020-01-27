#include "container/sv_db/query_objects/nucSeqSql.h"

using namespace libMA;

#ifdef WITH_PYTHON

void exportNucSeqSql( py::module& rxPyModuleId )
{
    // export the NucSeqFromSql classes
#ifndef USE_NEW_DB_API
    exportModule<AllNucSeqFromSql, std::shared_ptr<SV_DB>, int64_t, size_t, size_t>( rxPyModuleId, "AllNucSeqFromSql" );
    exportModule<NucSeqFromSql, std::shared_ptr<SV_DB>, int64_t>( rxPyModuleId, "NucSeqFromSql" );
    exportModule<PairedNucSeqFromSql, std::shared_ptr<SV_DB>, int64_t>( rxPyModuleId, "PairedNucSeqFromSql" );
#else
	using DBCon = MySQLConDB;
	exportModule<AllNucSeqFromSql<DBCon>, std::shared_ptr<_SV_DB<DBCon>>, int64_t, size_t, size_t>(rxPyModuleId, "AllNucSeqFromSql");
	exportModule<NucSeqFromSql<DBCon>, std::shared_ptr<_SV_DB<DBCon>>, int64_t>(rxPyModuleId, "NucSeqFromSql");
	exportModule<PairedNucSeqFromSql<DBCon>, std::shared_ptr<_SV_DB<DBCon>>, int64_t>(rxPyModuleId, "PairedNucSeqFromSql");
#endif
} // function

#endif