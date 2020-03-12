#include "msv/container/sv_db/query_objects/fetchSvJump.h"

using namespace libMSV;

#ifdef WITH_PYTHON

#include "ms/container/sv_db/py_db_conf.h"

void exportSvJump( libMS::SubmoduleOrganizer& xOrganizer )
{
    // export the SortedSvJumpFromSql class
    py::class_<SortedSvJumpFromSql<DBConSingle>>( xOrganizer.util( ), "SortedSvJumpFromSql" )
        .def( py::init<std::shared_ptr<DBConSingle>, int64_t>( ) )
        .def( "has_next_start", &SortedSvJumpFromSql<DBConSingle>::hasNextStart )
        .def( "has_next_end", &SortedSvJumpFromSql<DBConSingle>::hasNextEnd )
        .def( "next_start_is_smaller", &SortedSvJumpFromSql<DBConSingle>::nextStartIsSmaller )
        .def( "get_next_start", &SortedSvJumpFromSql<DBConSingle>::getNextStart )
        .def( "get_next_end", &SortedSvJumpFromSql<DBConSingle>::getNextEnd );

    py::class_<SvJumpFromSql<DBConSingle>>( xOrganizer.util( ), "SvJumpFromSql" )
        .def( py::init<std::shared_ptr<DBConSingle>, int64_t, int64_t, int64_t, uint32_t, uint32_t>( ) )
        .def( "has_next", &SvJumpFromSql<DBConSingle>::hasNext )
        .def( "get_next", &SvJumpFromSql<DBConSingle>::getNext );
} // function

#endif // WITH_PYTHON
