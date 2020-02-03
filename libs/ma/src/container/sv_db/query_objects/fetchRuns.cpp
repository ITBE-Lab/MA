#include "container/sv_db/query_objects/fetchRuns.h"

using namespace libMA;

#ifdef WITH_PYTHON

#include "container/sv_db/py_db_conf.h"

void exportRunsFromDb( py::module& rxPyModuleId )
{

    // export the SvCallerRunsFromDb class
    py::class_<SvCallerRunsFromDb<DBCon>>( rxPyModuleId, "SvCallerRunsFromDb" )
        .def( py::init<std::shared_ptr<SV_Schema<DBCon>>>( ) )
        .def( "id", &SvCallerRunsFromDb<DBCon>::id )
        .def( "name", &SvCallerRunsFromDb<DBCon>::name )
        .def( "desc", &SvCallerRunsFromDb<DBCon>::desc )
        .def( "eof", &SvCallerRunsFromDb<DBCon>::eof )
        .def( "next", &SvCallerRunsFromDb<DBCon>::next );
} // function
#endif
