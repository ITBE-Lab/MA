#include "container/sv_db/query_objects/callInserter.h"

using namespace libMA;

#ifdef WITH_PYTHON

void exportSvCallInserter( py::module& rxPyModuleId )
{
    // export the SvCallInserter class
    py::class_<SvCallInserter, std::shared_ptr<SvCallInserter>>( rxPyModuleId, "SvCallInserter" )
        .def( py::init<std::shared_ptr<SV_DB>, int64_t>( ) )
        .def( py::init<std::shared_ptr<SV_DB>, std::string, std::string, int64_t>( ) )
        .def_readonly( "sv_caller_run_id", &SvCallInserter::iSvCallerRunId )
        .def( "insert_call", &SvCallInserter::insertCall );
} // function

#endif