#include "ms/container/sv_db/pool_container.h"

using namespace libMS;

#ifdef WITH_PYTHON
#include "ms/container/sv_db/py_db_conf.h"


void exportPoolContainer( py::module& rxPyModuleId )
{
    py::class_<PoolContainer<DBCon>, Container, std::shared_ptr<PoolContainer<DBCon>>>( rxPyModuleId,
                                                                                                  "PoolContainer" )
        .def( py::init<size_t, std::string>( ) );

} // function
#endif // WITH_PYTHON
