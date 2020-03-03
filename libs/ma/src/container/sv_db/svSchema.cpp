#include "container/sv_db/svSchema.h"
#include "container/container.h"
#include "container/sv_db/py_db_conf.h"
#include "module/combineOverlappingCalls.h"

// include classes that implement sql queries
#include "container/sv_db/query_objects/callInserter.h"
#include "container/sv_db/query_objects/fetchCalls.h"
#include "container/sv_db/query_objects/fetchSvJump.h"
#include "container/sv_db/query_objects/jumpInserter.h"
#include "container/sv_db/query_objects/nucSeqSql.h"
#include "container/sv_db/query_objects/readInserter.h"

using namespace libMA;


uint32_t getNumJumpsInArea( std::shared_ptr<DBConSingle> pConnection, std::shared_ptr<Pack> pPack, int64_t iRunId,
                            int64_t iX, int64_t iY, uint64_t uiW, uint64_t uiH, uint64_t uiLimit )
{
    uint32_t uiX = 0;
    if( iX > 0 )
        uiX = (uint32_t)iX;
    uint32_t uiY = 0;
    if( iY > 0 )
        uiY = (uint32_t)iY;
    if( uiX + uiW > pPack->uiUnpackedSizeForwardStrand )
        uiW = pPack->uiUnpackedSizeForwardStrand - uiX;
    if( uiY + uiH > pPack->uiUnpackedSizeForwardStrand )
        uiH = pPack->uiUnpackedSizeForwardStrand - uiY;

    SQLQuery<DBConSingle, uint32_t> xQuery( pConnection,
                                            "SELECT COUNT(*) "
                                            "FROM sv_jump_table "
                                            "WHERE sv_jump_run_id = ? "
                                            "AND ( (from_pos >= ? AND from_pos <= ?) OR from_pos = ? ) "
                                            "AND ( (to_pos >= ? AND to_pos <= ?) OR to_pos = ? ) "
                                            "LIMIT ? " );
    // FIXME: Don't use scalar anymore!
    return xQuery.scalar( iRunId, uiX, uiX + (uint32_t)uiW, std::numeric_limits<uint32_t>::max( ), uiY,
                          uiY + (uint32_t)uiH, std::numeric_limits<uint32_t>::max( ), uiLimit );
} // function


#ifdef WITH_PYTHON
#include "pybind11_json/pybind11_json.hpp"

void exportSoCDbWriter( py::module& rxPyModuleId )
{
    py::class_<DBConSingle, std::shared_ptr<DBConSingle>>( rxPyModuleId, "DbConn" )
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

    py::class_<SvCallTable<DBConSingle>, std::shared_ptr<SvCallTable<DBConSingle>>>( rxPyModuleId, "SvCallTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "reconstruct_sequenced_genome", &SvCallTable<DBConSingle>::reconstructSequencedGenome )
        .def( "num_calls", &SvCallTable<DBConSingle>::numCalls_py )
        .def( "max_score", &SvCallTable<DBConSingle>::maxScore )
        .def( "min_score", &SvCallTable<DBConSingle>::minScore )
        .def( "call_area", &SvCallTable<DBConSingle>::callArea )
        .def( "drop_indices", &SvCallTable<DBConSingle>::dropIndices )
        .def( "gen_indices", &SvCallTable<DBConSingle>::genIndices );

    using X = SvCallTableAnalyzer<DBCon, false>;
    py::class_<X, std::shared_ptr<X>>( rxPyModuleId, "SvCallTableAnalyzer" )
        .def( py::init<std::shared_ptr<PoolContainer<DBCon>>>( ) )
        .def( "num_overlaps", &X::numOverlaps )
        .def( "num_invalid_calls", &X::numInvalidCalls )
        .def( "blur_on_overlaps", &X::blurOnOverlaps );

    py::class_<SvJumpRunTable<DBConSingle>, std::shared_ptr<SvJumpRunTable<DBConSingle>>>( rxPyModuleId,
                                                                                           "JumpRunTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) );

    py::class_<SvJumpTable<DBConSingle>, std::shared_ptr<SvJumpTable<DBConSingle>>>( rxPyModuleId, "SvJumpTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "create_indices", &SvJumpTable<DBConSingle>::createIndices )
        .def( "drop_indices", &SvJumpTable<DBConSingle>::dropIndices )
        .def( "num_jumps", &SvJumpTable<DBConSingle>::numJumps );

    py::class_<ReadTable<DBConSingle>, std::shared_ptr<ReadTable<DBConSingle>>>( rxPyModuleId, "ReadTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "get_read", &ReadTable<DBConSingle>::getRead );

    py::class_<SvCallSupportTable<DBConSingle>, std::shared_ptr<SvCallSupportTable<DBConSingle>>>(
        rxPyModuleId, "SvCallSupportTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) );

    py::class_<SvCallerRunTable<DBConSingle>, std::shared_ptr<SvCallerRunTable<DBConSingle>>>( rxPyModuleId,
                                                                                               "SvCallerRunTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "getIds", &SvCallerRunTable<DBConSingle>::getIds )
        .def( "getName", &SvCallerRunTable<DBConSingle>::getName )
        .def( "getDesc", &SvCallerRunTable<DBConSingle>::getDesc )
        .def( "jump_run_id", &SvCallerRunTable<DBConSingle>::getSvJumpRunId )
        .def( "newest_unique_runs", &SvCallerRunTable<DBConSingle>::getNewestUnique )
        .def( "getDate", &SvCallerRunTable<DBConSingle>::getDate );

    py::class_<rect>( rxPyModuleId, "rect" ) //
        .def_readonly( "x", &rect::x )
        .def_readonly( "y", &rect::y )
        .def_readonly( "w", &rect::w )
        .def_readonly( "h", &rect::h )
        .def_readonly( "i", &rect::i )
        .def_readonly( "j", &rect::j )
        .def_readonly( "c", &rect::c );
    py::bind_vector<std::vector<rect>>( rxPyModuleId, "rectVector", "docstr" );
    rxPyModuleId.def( "get_num_jumps_in_area", &getNumJumpsInArea );
    rxPyModuleId.def( "get_call_overview", &getCallOverview<DBCon> );
    rxPyModuleId.def( "get_call_overview_area", &getCallOverviewArea<DBConSingle> );

    rxPyModuleId.def( "combine_overlapping_calls", &combineOverlappingCalls<DBCon> );

    exportSvCallInserter( rxPyModuleId );
    exportCallsFromDb( rxPyModuleId );
    exportSvJump( rxPyModuleId );
    exportSvJumpInserter( rxPyModuleId );
    exportNucSeqSql( rxPyModuleId );
    exportReadInserter( rxPyModuleId );
} // function
#endif // WITH_PYTHON
