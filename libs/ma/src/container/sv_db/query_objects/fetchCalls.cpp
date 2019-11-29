#include "container/sv_db/query_objects/fetchCalls.h"

using namespace libMA;

#ifdef WITH_PYTHON

void exportCallsFromDb( py::module& rxPyModuleId )
{
    // export the SvCallsFromDb class
    py::class_<SvCallsFromDb>( rxPyModuleId, "SvCallsFromDb" )
        .def( py::init<const ParameterSetManager&, std::shared_ptr<SV_DB>, int64_t>( ) )
        .def( py::init<const ParameterSetManager&, std::shared_ptr<SV_DB>, int64_t, double>( ) )
        .def( py::init<const ParameterSetManager&, std::shared_ptr<SV_DB>, int64_t, int64_t, bool, int64_t>( ) )
        .def( py::init<const ParameterSetManager&, std::shared_ptr<SV_DB>, int64_t, uint32_t, uint32_t, uint32_t,
                       uint32_t>( ) )
        .def( py::init<const ParameterSetManager&, std::shared_ptr<SV_DB>, int64_t, int64_t, int64_t, int64_t, int64_t,
                       double>( ) )
        .def( "next", &SvCallsFromDb::next )
        .def( "hasNext", &SvCallsFromDb::hasNext );
} // function

#endif