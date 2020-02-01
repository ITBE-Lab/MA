#include "container/sv_db/query_objects/nucSeqSql.h"

using namespace libMA;

#ifdef WITH_PYTHON

void exportNucSeqSql( py::module& rxPyModuleId )
{
    // export the NucSeqFromSql classes
    using DBCon = SQLDB<MySQLConDB>;
	exportModule<AllNucSeqFromSql<DBCon>, std::shared_ptr<_SV_DB<DBCon>>, int64_t, size_t, size_t>(rxPyModuleId, "AllNucSeqFromSql");
	exportModule<NucSeqFromSql<DBCon>, std::shared_ptr<_SV_DB<DBCon>>, int64_t>(rxPyModuleId, "NucSeqFromSql");
	exportModule<PairedNucSeqFromSql<DBCon>, std::shared_ptr<_SV_DB<DBCon>>, int64_t>(rxPyModuleId, "PairedNucSeqFromSql");
} // function

#endif