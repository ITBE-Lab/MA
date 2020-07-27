#include "msv/container/sv_db/query_objects/kMerInserter.h"

using namespace libMSV;

#ifdef WITH_PYTHON

#include "ms/container/sv_db/py_db_conf.h"


void exportKMerInserter( libMS::SubmoduleOrganizer& xOrganizer )
{
    exportInserterContainer<GetKMerInserterContainerModule<DBCon, DBConSingle>>( xOrganizer, "KMerInserter" );

    exportModule<KMerInserterModule<DBCon>>( xOrganizer, "KMerInserterModule" );
} // function

#endif