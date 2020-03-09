#include "msv/container/sv_db/query_objects/callInserter.h"

using namespace libMSV;

#ifdef WITH_PYTHON

#include "ms/container/sv_db/py_db_conf.h"

void exportSvCallInserter( libMS::SubmoduleOrganizer& xOrganizer )
{
    // export the CallInserter class
    exportInserterContainer<GetCallInserterContainerModule<DBCon, DBConSingle>>( xOrganizer, "CallInserter" );

    exportInserterContainer<GetCallVectorInserterContainerModule<DBCon, DBConSingle>>( xOrganizer,
                                                                                       "CallVectorInserter" );

    exportModule<SvCallInserterModule<DBCon>>( xOrganizer, "CallInserterModule" );
    exportModule<SvCallVectorInserterModule<DBCon>>( xOrganizer, "CallVectorInserterModule" );
} // function

#endif // WITH_PYTHON
