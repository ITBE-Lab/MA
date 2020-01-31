#include "container/sv_db/query_objects/fetchRuns.h"

using namespace libMA;

#ifdef WITH_PYTHON

#ifndef USE_NEW_DB_API
void exportRunsFromDb( py::module& rxPyModuleId )
{

    // export the SvCallerRunsFromDb class
    py::class_<SvCallerRunsFromDb>( rxPyModuleId, "SvCallerRunsFromDb" )
        .def( py::init<std::shared_ptr<SV_DB>>( ) )
        .def( "id", &SvCallerRunsFromDb::id )
        .def( "name", &SvCallerRunsFromDb::name )
        .def( "desc", &SvCallerRunsFromDb::desc )
        .def( "eof", &SvCallerRunsFromDb::eof )
        .def( "next", &SvCallerRunsFromDb::next );
} // function
#else

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
#endif // USE_NEW_DB_API
#endif
