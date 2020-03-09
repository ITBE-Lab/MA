#include "msv/container/sv_db/query_objects/jumpInserter.h"

using namespace libMSV;

#ifdef WITH_PYTHON

#include "ms/container/sv_db/py_db_conf.h"

void exportSvJumpInserter( libMS::SubmoduleOrganizer& xOrganizer )
{
    exportInserterContainer<GetJumpInserterContainerModule<DBCon, DBConSingle>>( xOrganizer, "JumpInserter" );

    exportModule<JumpInserterModule<DBCon>>( xOrganizer, "JumpInserterModule" );
} // function

#endif
