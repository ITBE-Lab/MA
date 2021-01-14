#include "ms/container/sv_db/pool_container.h"


#ifdef WITH_DB

using namespace libMS;

#ifdef WITH_PYTHON
#include "ms/container/sv_db/py_db_conf.h"
#include "pybind11_json/pybind11_json.hpp"

void exportPoolContainer( SubmoduleOrganizer& xOrganizer )
{
    py::class_<PoolContainer<DBCon>, Container, std::shared_ptr<PoolContainer<DBCon>>>( xOrganizer.container( ),
                                                                                        "PoolContainer" )
        .def( py::init<size_t, std::string>( ) );

    py::class_<DBConSingle, std::shared_ptr<DBConSingle>>( xOrganizer.util( ), "DbConn" )
        .def( py::init<std::string>( ) )
        /* This makes it so, that DbConn can be initialized from a python dictionary.
         * It makes use of https://github.com/pybind/pybind11_json
         * For some reason the nlohmann::json object can not be passed directly to py::init,
         * however the py::object is converted automatically since the header pybind11_json.hpp is included here.
         * @todo we could drop the guy an issue asking/suggesting to make it possible to putt in the json directly,
         * which would make the code more readable
         */
        .def( py::init<py::object /* = json */>( ) )
        .def( "drop_schema", &DBConSingle::dropSchema );
} // function
#endif // WITH_PYTHON

#endif // WITH_DB
