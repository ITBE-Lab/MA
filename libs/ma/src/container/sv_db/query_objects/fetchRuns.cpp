#include "container/sv_db/query_objects/fetchRuns.h"

using namespace libMA;

#ifdef WITH_PYTHON

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

#endif