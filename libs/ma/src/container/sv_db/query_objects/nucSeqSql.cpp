#include "container/sv_db/query_objects/nucSeqSql.h"

using namespace libMA;

#ifdef WITH_PYTHON

void exportNucSeqSql( py::module& rxPyModuleId )
{
    // export the NucSeqFromSql classes
    exportModule<AllNucSeqFromSql, std::shared_ptr<SV_DB>, int64_t, size_t, size_t>( rxPyModuleId, "AllNucSeqFromSql" );
    exportModule<NucSeqFromSql, std::shared_ptr<SV_DB>, int64_t>( rxPyModuleId, "NucSeqFromSql" );
    exportModule<PairedNucSeqFromSql, std::shared_ptr<SV_DB>, int64_t>( rxPyModuleId, "PairedNucSeqFromSql" );
} // function

#endif