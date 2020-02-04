#include "container/sv_db/query_objects/nucSeqSql.h"

using namespace libMA;

#ifdef WITH_PYTHON

#include "container/sv_db/py_db_conf.h"

void exportNucSeqSql( py::module& rxPyModuleId )
{
#if 0 // @todo
    // export the NucSeqFromSql classes
	exportModule<AllNucSeqFromSql<DBCon>, std::shared_ptr<SV_Schema<DBCon>>, int64_t, size_t, size_t>(rxPyModuleId, "AllNucSeqFromSql");
	exportModule<NucSeqFromSql<DBCon>, std::shared_ptr<SV_Schema<DBCon>>, int64_t>(rxPyModuleId, "NucSeqFromSql");
	exportModule<PairedNucSeqFromSql<DBCon>, std::shared_ptr<SV_Schema<DBCon>>, int64_t>(rxPyModuleId, "PairedNucSeqFromSql");
#endif
} // function

#endif