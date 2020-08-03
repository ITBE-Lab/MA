#include "msv/container/sv_db/query_objects/readInserter.h"

using namespace libMSV;

#ifdef WITH_PYTHON

#include "ms/container/sv_db/py_db_conf.h"


void exportReadInserter( libMS::SubmoduleOrganizer& xOrganizer )
{
    exportInserterContainer<GetReadInserterContainerModule<DBCon, DBConSingle>>( xOrganizer, "ReadInserter" );
    exportInserterContainer<GetPairedReadInserterContainerModule<DBCon, DBConSingle>>( xOrganizer, "PairedReadInserter" );

    exportModule<ReadInserterModule<DBCon>>( xOrganizer, "ReadInserterModule" );
    exportModule<PairedReadInserterModule<DBCon>>( xOrganizer, "PairedReadInserterModule" );
} // function

#endif