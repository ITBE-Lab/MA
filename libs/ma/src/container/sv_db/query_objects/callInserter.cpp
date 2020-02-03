#include "container/sv_db/query_objects/callInserter.h"

using namespace libMA;

#ifdef WITH_PYTHON

#include "container/sv_db/py_db_conf.h"

void exportSvCallInserter( py::module& rxPyModuleId )
{
    // export the SvCallInserter class
    py::class_<SvCallInserter<DBCon>, std::shared_ptr<SvCallInserter<DBCon>>>( rxPyModuleId, "SvCallInserter" )
        .def( py::init<std::shared_ptr<DBCon>, int64_t>( ) )
        .def( py::init<std::shared_ptr<DBCon>, std::string, std::string, int64_t>( ) )
        .def_readonly( "sv_caller_run_id", &SvCallInserter<DBCon>::iSvCallerRunId )
        .def( "insert_call", &SvCallInserter<DBCon>::insertCall )
        .def( "end_transaction", &SvCallInserter<DBCon>::endTransaction );
} // function

#endif // WITH_PYTHON