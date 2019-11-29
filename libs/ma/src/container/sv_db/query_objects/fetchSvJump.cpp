#include "container/sv_db/query_objects/fetchSvJump.h"

using namespace libMA;

#ifdef WITH_PYTHON

void exportSvJump( py::module& rxPyModuleId )
{
    // export the SortedSvJumpFromSql class
    py::class_<SortedSvJumpFromSql>( rxPyModuleId, "SortedSvJumpFromSql" )
        .def( py::init<const ParameterSetManager&, std::shared_ptr<SV_DB>, int64_t>( ) )
        .def( py::init<const ParameterSetManager&, std::shared_ptr<SV_DB>, int64_t, int64_t, int64_t, uint32_t,
                       uint32_t>( ) )
        .def( "has_next_start", &SortedSvJumpFromSql::hasNextStart )
        .def( "has_next_end", &SortedSvJumpFromSql::hasNextEnd )
        .def( "next_start_is_smaller", &SortedSvJumpFromSql::nextStartIsSmaller )
        .def( "get_next_start", &SortedSvJumpFromSql::getNextStart )
        .def( "get_next_end", &SortedSvJumpFromSql::getNextEnd );
} // function

#endif