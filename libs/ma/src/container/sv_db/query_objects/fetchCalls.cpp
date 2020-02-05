#include "container/sv_db/query_objects/fetchCalls.h"

using namespace libMA;

#ifdef WITH_PYTHON

#include "container/sv_db/py_db_conf.h"
void exportCallsFromDb( py::module& rxPyModuleId )
{
    // export the SvCallsFromDb class
    py::class_<SvCallsFromDb<DBCon>>( rxPyModuleId, "SvCallsFromDb" )
        .def( py::init<const ParameterSetManager&, std::shared_ptr<DBCon>, int64_t>( ) )
        .def( py::init<const ParameterSetManager&, std::shared_ptr<DBCon>, int64_t, double>( ) )
        .def( py::init<const ParameterSetManager&, std::shared_ptr<DBCon>, int64_t, int64_t, bool, int64_t>( ) )
        .def( py::init<const ParameterSetManager&, std::shared_ptr<DBCon>, int64_t, uint32_t, uint32_t, uint32_t,
                       uint32_t>( ) )
        .def( py::init<const ParameterSetManager&, std::shared_ptr<DBCon>, int64_t, int64_t, int64_t, int64_t, int64_t,
                       double>( ) )
        .def( "next", &SvCallsFromDb<DBCon>::next )
        .def( "hasNext", &SvCallsFromDb<DBCon>::hasNext );
} // function

#endif
