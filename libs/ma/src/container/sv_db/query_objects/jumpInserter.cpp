#include "container/sv_db/query_objects/jumpInserter.h"

using namespace libMA;

#ifdef WITH_PYTHON

void exportSvJumpInserter( py::module& rxPyModuleId )
{
    // export the SvJumpInserter class
    py::class_<SvJumpInserter, std::shared_ptr<SvJumpInserter>>( rxPyModuleId, "SvJumpInserter" )
        .def( py::init<std::shared_ptr<SV_DB>, std::string, std::string>( ) )
        .def( "read_context", &SvJumpInserter::readContext )
        .def_readonly( "sv_jump_run_id", &SvJumpInserter::iSvJumpRunId );

    // export the ReadContex class
    py::class_<SvJumpInserter::ReadContex>( rxPyModuleId, "ReadContex" )
        .def( "insert_jump", &SvJumpInserter::ReadContex::insertJump );
    exportModule<SvDbInserter, std::shared_ptr<SV_DB>, std::string>(
        rxPyModuleId, "SvDbInserter", []( auto&& x ) { x.def_readonly( "jump_inserter", &SvDbInserter::xInserter ); } );
    exportModule<BufferedSvDbInserter, std::shared_ptr<SV_DB>, int64_t>(
        rxPyModuleId, "BufferedSvDbInserter", []( auto&& x ) { x.def( "commit", &BufferedSvDbInserter::commit ); } );
} // function

#endif