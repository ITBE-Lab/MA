#include "msv/container/sv_db/query_objects/fetchCalls.h"

using namespace libMSV;

#ifdef WITH_PYTHON

#include <pybind11/stl.h>
#include "ms/container/sv_db/py_db_conf.h"
void exportCallsFromDb( libMS::SubmoduleOrganizer& xOrganizer )
{
    // export the SvCallsFromDb class
    py::class_<SvCallsFromDb<DBConSingle>>( xOrganizer.util( ), "SvCallsFromDb" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "init", ( void ( SvCallsFromDb<DBConSingle>::* )( int64_t, int64_t, int64_t, int64_t, int64_t, int64_t,
                                                                bool, int64_t, double, double ) ) &
                          SvCallsFromDb<DBConSingle>::initFetchQuery )
        .def( "init", ( void ( SvCallsFromDb<DBConSingle>::* )( double, double, int64_t, int64_t, int64_t, int64_t,
                                                                int64_t, int64_t, bool, int64_t ) ) &
                          SvCallsFromDb<DBConSingle>::initFetchQuery )
        .def( "init", ( void ( SvCallsFromDb<DBConSingle>::* )( int64_t, int64_t, int64_t, int64_t, int64_t, int64_t,
                                                                bool, int64_t ) ) &
                          SvCallsFromDb<DBConSingle>::initFetchQuery )
        .def( "init", ( void ( SvCallsFromDb<DBConSingle>::* )( int64_t, int64_t, int64_t, int64_t, int64_t, double,
                                                                double ) ) &
                          SvCallsFromDb<DBConSingle>::initFetchQuery )
        .def( "init", ( void ( SvCallsFromDb<DBConSingle>::* )( int64_t, int64_t, int64_t, int64_t, int64_t ) ) &
                          SvCallsFromDb<DBConSingle>::initFetchQuery )
        .def( "next", &SvCallsFromDb<DBConSingle>::next )
        .def( "count", &SvCallsFromDb<DBConSingle>::count )
        .def( "hasNext", &SvCallsFromDb<DBConSingle>::hasNext );
} // function

#endif
