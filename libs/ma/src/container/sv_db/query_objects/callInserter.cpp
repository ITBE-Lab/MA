#include "container/sv_db/query_objects/callInserter.h"

using namespace libMA;

#ifdef WITH_PYTHON
#ifndef USE_NEW_DB_API
void exportSvCallInserter( py::module& rxPyModuleId )
{
    // export the SvCallInserter class
    py::class_<SvCallInserter, std::shared_ptr<SvCallInserter>>( rxPyModuleId, "SvCallInserter" )
        .def( py::init<std::shared_ptr<SV_DB>, int64_t>( ) )
        .def( py::init<std::shared_ptr<SV_DB>, std::string, std::string, int64_t>( ) )
        .def_readonly( "sv_caller_run_id", &SvCallInserter::iSvCallerRunId )
        .def( "insert_call", &SvCallInserter::insertCall )
        .def( "end_transaction", &SvCallInserter::endTransaction );
} // function
#else

/* NEW DATABASE INTERFACE */

using DBCon = MySQLConDB;

void exportSvCallInserter( py::module& rxPyModuleId )
{
    // export the SvCallInserter class
    py::class_<SvCallInserter<DBCon>, std::shared_ptr<SvCallInserter<DBCon>>>( rxPyModuleId, "SvCallInserter" )
        .def( py::init<std::shared_ptr<_SV_DB<DBCon>>, int64_t>( ) )
        .def( py::init<std::shared_ptr<_SV_DB<DBCon>>, std::string, std::string, int64_t>( ) )
        .def_readonly( "sv_caller_run_id", &SvCallInserter<DBCon>::iSvCallerRunId )
        .def( "insert_call", &SvCallInserter<DBCon>::insertCall )
        .def( "end_transaction", &SvCallInserter<DBCon>::endTransaction );
} // function
#endif // USE_NEW_DB_API

#endif // WITH_PYTHON