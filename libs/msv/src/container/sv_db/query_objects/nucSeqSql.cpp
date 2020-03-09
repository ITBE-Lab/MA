#include "msv/container/sv_db/query_objects/nucSeqSql.h"

#ifdef WITH_PYTHON

#include "ms/container/sv_db/py_db_conf.h"

using namespace libMSV;
using namespace libMS;

void exportNucSeqSql( libMS::SubmoduleOrganizer& xOrganizer )
{
    py::class_<NucSeqQueryContainer<DBCon>, Container, std::shared_ptr<NucSeqQueryContainer<DBCon>>>(
        xOrganizer.container( ), "NucSeqQueryContainer" );

    exportModule<GetNucSeqFromSqlQuery<DBCon>, int64_t, size_t, size_t, bool, bool>( xOrganizer,
                                                                                     "GetNucSeqFromSqlQuery" );

    exportModule<NucSeqFetcher<DBCon>>( xOrganizer, "NucSeqFetcher" );
} // function

#endif