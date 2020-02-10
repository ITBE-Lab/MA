#include "container/sv_db/query_objects/callInserter.h"

using namespace libMA;

#ifdef WITH_PYTHON

#include "container/sv_db/py_db_conf.h"

void exportSvCallInserter( py::module& rxPyModuleId )
{
    // export the CallInserter class
    exportInserterContainer<GetCallInserterContainerModule<DBCon, DBConSingle>>
        ( rxPyModuleId, "CallInserter" );

    exportInserterContainer<GetCallVectorInserterContainerModule<DBCon, DBConSingle>>
        ( rxPyModuleId, "CallVectorInserter" );

    exportModule<SvCallInserterModule<DBCon>>( rxPyModuleId, "CallInserterModule" );
    exportModule<SvCallVectorInserterModule<DBCon>>( rxPyModuleId, "CallVectorInserterModule" );
} // function

#endif // WITH_PYTHON
