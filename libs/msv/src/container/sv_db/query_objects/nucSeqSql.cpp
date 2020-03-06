#include "msv/container/sv_db/query_objects/nucSeqSql.h"

#ifdef WITH_PYTHON

#include "ms/container/sv_db/py_db_conf.h"

using namespace libMSV;
using namespace libMS;

void exportNucSeqSql( py::module& rxPyModuleId )
{
    py::class_<NucSeqQueryContainer<DBCon>, Container, std::shared_ptr<NucSeqQueryContainer<DBCon>>>(
        rxPyModuleId, "NucSeqQueryContainer" );

    exportModule<GetNucSeqFromSqlQuery<DBCon>, int64_t, size_t, size_t, bool, bool>( rxPyModuleId,
                                                                                     "GetNucSeqFromSqlQuery" );

    exportModule<NucSeqFetcher<DBCon>>( rxPyModuleId, "NucSeqFetcher" );
} // function

#endif