#include "container/sv_db/query_objects/jumpInserter.h"

using namespace libMA;

#ifdef WITH_PYTHON
#ifndef USE_NEW_DB_API
void exportSvJumpInserter( py::module& rxPyModuleId )
{
    // export the SvJumpInserter class
    py::class_<SvJumpInserter, std::shared_ptr<SvJumpInserter>>( rxPyModuleId, "SvJumpInserter" )
        .def( py::init<std::shared_ptr<SV_DB>, std::string, std::string>( ) )
        .def( "read_context", &SvJumpInserter::readContext )
        .def( "end_transaction", &SvJumpInserter::endTransaction )
        .def_readonly( "sv_jump_run_id", &SvJumpInserter::iSvJumpRunId );

    // export the ReadContex class
    py::class_<SvJumpInserter::ReadContex>( rxPyModuleId, "ReadContex" )
        .def( "insert_jump", &SvJumpInserter::ReadContex::insertJump );
    exportModule<SvDbInserter, std::shared_ptr<SV_DB>, std::string>(
        rxPyModuleId, "SvDbInserter", []( auto&& x ) { x.def_readonly( "jump_inserter", &SvDbInserter::xInserter ); } );
    exportModule<BufferedSvDbInserter, std::shared_ptr<SvJumpInserter>>(
        rxPyModuleId, "BufferedSvDbInserter", []( auto&& x ) { x.def( "commit", &BufferedSvDbInserter::commit ); } );
} // function
#else

using DBCon = MySQLConDB;

void exportSvJumpInserter( py::module& rxPyModuleId )
{
    // export the SvJumpInserter class
    py::class_<SvJumpInserter<DBCon>, std::shared_ptr<SvJumpInserter<DBCon>>>( rxPyModuleId, "SvJumpInserter" )
        .def( py::init<std::shared_ptr<_SV_DB<DBCon>>, std::string, std::string>( ) )
        .def( "read_context", &SvJumpInserter<DBCon>::readContext )
        .def( "end_transaction", &SvJumpInserter<DBCon>::endTransaction )
        .def_readonly( "sv_jump_run_id", &SvJumpInserter<DBCon>::iSvJumpRunId );

    // export the ReadContex class
    py::class_<SvJumpInserter<DBCon>::ReadContex>( rxPyModuleId, "ReadContex" )
        .def( "insert_jump", &SvJumpInserter<DBCon>::ReadContex::insertJump );
    exportModule<SvDbInserter<DBCon>, std::shared_ptr<_SV_DB<DBCon>>, std::string>(
        rxPyModuleId, "SvDbInserter",
        []( auto&& x ) { x.def_readonly( "jump_inserter", &SvDbInserter<DBCon>::xInserter ); } );
    exportModule<BufferedSvDbInserter<DBCon>, std::shared_ptr<SvJumpInserter<DBCon>>>(
        rxPyModuleId, "BufferedSvDbInserter",
        []( auto&& x ) { x.def( "commit", &BufferedSvDbInserter<DBCon>::commit ); } );
} // function

#endif // USE_NEW_DB_API
#endif
