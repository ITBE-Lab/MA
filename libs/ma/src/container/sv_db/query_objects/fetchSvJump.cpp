#include "container/sv_db/query_objects/fetchSvJump.h"

using namespace libMA;

#ifdef WITH_PYTHON

#include "container/sv_db/py_db_conf.h"

void exportSvJump( py::module& rxPyModuleId )
{
    // export the SortedSvJumpFromSql class
    py::class_<SortedSvJumpFromSql<DBConSingle>>( rxPyModuleId, "SortedSvJumpFromSql" )
        .def( py::init<std::shared_ptr<DBConSingle>, int64_t>( ) )
        .def( py::init<std::shared_ptr<DBConSingle>, int64_t, int64_t, int64_t, uint32_t,
                       uint32_t>( ) )
        .def( "has_next_start", &SortedSvJumpFromSql<DBConSingle>::hasNextStart )
        .def( "has_next_end", &SortedSvJumpFromSql<DBConSingle>::hasNextEnd )
        .def( "next_start_is_smaller", &SortedSvJumpFromSql<DBConSingle>::nextStartIsSmaller )
        .def( "get_next_start", &SortedSvJumpFromSql<DBConSingle>::getNextStart )
        .def( "get_next_end", &SortedSvJumpFromSql<DBConSingle>::getNextEnd );
} // function

#endif // WITH_PYTHON
