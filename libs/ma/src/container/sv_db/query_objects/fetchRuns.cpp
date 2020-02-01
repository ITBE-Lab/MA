#include "container/sv_db/query_objects/fetchRuns.h"

using namespace libMA;

#ifdef WITH_PYTHON

using DBCon = SQLDB<MySQLConDB>;

void exportRunsFromDb( py::module& rxPyModuleId )
{

    // export the SvCallerRunsFromDb class
    py::class_<SvCallerRunsFromDb<DBCon>>( rxPyModuleId, "SvCallerRunsFromDb" )
        .def( py::init<std::shared_ptr<_SV_DB<DBCon>>>( ) )
        .def( "id", &SvCallerRunsFromDb<DBCon>::id )
        .def( "name", &SvCallerRunsFromDb<DBCon>::name )
        .def( "desc", &SvCallerRunsFromDb<DBCon>::desc )
        .def( "eof", &SvCallerRunsFromDb<DBCon>::eof )
        .def( "next", &SvCallerRunsFromDb<DBCon>::next );
} // function
#endif
