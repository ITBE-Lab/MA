#include "container/sv_db/query_objects/callInserter.h"

using namespace libMA;

#ifdef WITH_PYTHON

#include "container/sv_db/py_db_conf.h"

void exportSvCallInserter( py::module& rxPyModuleId )
{
    // export the CallInserter class
    exportInserterContainer<GetCallInserterContainerModule<DBCon, DBConSingle>,
                            std::string, // name
                            std::string, // desc
                            int64_t, // timestamp
                            int64_t> // sv_jump_run_id
        ( rxPyModuleId, "CallInserter" );

    exportModule<CallInserterModule>( rxPyModuleId, "CallInserterModule" );
} // function

#endif // WITH_PYTHON