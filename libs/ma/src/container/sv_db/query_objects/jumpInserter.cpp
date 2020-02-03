#include "container/sv_db/query_objects/jumpInserter.h"

using namespace libMA;

#ifdef WITH_PYTHON

#include "container/sv_db/py_db_conf.h"

void exportSvJumpInserter( py::module& rxPyModuleId )
{
    // export the SvJumpInserter class
    py::class_<SvJumpInserter<DBCon>, std::shared_ptr<SvJumpInserter<DBCon>>>( rxPyModuleId, "SvJumpInserter" )
        .def( py::init<std::shared_ptr<SV_Schema<DBCon>>, std::string, std::string>( ) )
        .def( "read_context", &SvJumpInserter<DBCon>::readContext )
        .def( "end_transaction", &SvJumpInserter<DBCon>::endTransaction )
        .def_readonly( "sv_jump_run_id", &SvJumpInserter<DBCon>::iSvJumpRunId );

    // export the ReadContex class
    py::class_<SvJumpInserter<DBCon>::ReadContex>( rxPyModuleId, "ReadContex" )
        .def( "insert_jump", &SvJumpInserter<DBCon>::ReadContex::insertJump );
    exportModule<SvDbInserter<DBCon>, std::shared_ptr<SV_Schema<DBCon>>, std::string>(
        rxPyModuleId, "SvDbInserter",
        []( auto&& x ) { x.def_readonly( "jump_inserter", &SvDbInserter<DBCon>::xInserter ); } );
    exportModule<BufferedSvDbInserter<DBCon>, std::shared_ptr<SvJumpInserter<DBCon>>>(
        rxPyModuleId, "BufferedSvDbInserter",
        []( auto&& x ) { x.def( "commit", &BufferedSvDbInserter<DBCon>::commit ); } );
} // function

#endif
