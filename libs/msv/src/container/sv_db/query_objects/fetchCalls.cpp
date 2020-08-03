#include "msv/container/sv_db/query_objects/fetchCalls.h"

using namespace libMSV;

#ifdef WITH_PYTHON

#include "ms/container/sv_db/py_db_conf.h"
void exportCallsFromDb( libMS::SubmoduleOrganizer& xOrganizer )
{
    // export the SvCallsFromDb class
    py::class_<SvCallsFromDb<DBConSingle>>( xOrganizer.util(), "SvCallsFromDb" )
        .def( py::init<const ParameterSetManager&, std::shared_ptr<DBConSingle>, int64_t>( ) )
        .def( py::init<const ParameterSetManager&, std::shared_ptr<DBConSingle>, int64_t, double>( ) )
        .def( py::init<const ParameterSetManager&, std::shared_ptr<DBConSingle>, int64_t, int64_t, bool, int64_t>( ) )
        .def( py::init<const ParameterSetManager&, std::shared_ptr<DBConSingle>, int64_t, uint32_t, uint32_t, uint32_t,
                       uint32_t>( ) )
        .def( py::init<const ParameterSetManager&, std::shared_ptr<DBConSingle>, int64_t, int64_t, int64_t, int64_t,
                       int64_t, double>( ) )
        .def( "next", &SvCallsFromDb<DBConSingle>::next )
        .def( "hasNext", &SvCallsFromDb<DBConSingle>::hasNext );
} // function

#endif
