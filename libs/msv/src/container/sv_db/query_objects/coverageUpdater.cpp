#include "msv/container/sv_db/query_objects/coverageUpdater.h"

using namespace libMSV;

#ifdef WITH_PYTHON

#include "ms/container/sv_db/py_db_conf.h"

void exportCoverageUpdater( libMS::SubmoduleOrganizer& xOrganizer )
{
    exportInserterContainer<GetCoverageUpdaterContainerModule<DBCon, DBConSingle>>( xOrganizer, "CoverageUpdater" );

    exportModule<CoverageUpdaterModule<DBCon>>( xOrganizer, "CoverageUpdaterModule" );
} // function

#endif
