#include "container/sv_db/svSchema.h"
#include "container/container.h"
#include "container/pack.h"
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

uint32_t getCallOverviewArea( std::shared_ptr<DBConSingle> pConnection, std::shared_ptr<Pack> pPack, int64_t iRunId,
                              double dMinScore, int64_t iX, int64_t iY, uint64_t uiW, uint64_t uiH )
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
                                            "FROM sv_call_table "
                                            "WHERE sv_caller_run_id = ? " // dim 1
                                            "AND ST_Overlaps(rectangle, ST_PolyFromWKB(?, 0)) "
                                            "AND score >= ? " );


    auto xWkb = WKBUint64Rectangle( geom::Rectangle<nucSeqIndex>( uiX, uiY, uiW, uiH ) );
    return xQuery.scalar( iRunId, xWkb, dMinScore );
} // function

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

struct rect
{
    uint32_t x, y, w, h, c, i, j;
    rect( uint32_t x, uint32_t y, uint32_t w, uint32_t h, uint32_t c, uint32_t i, uint32_t j )
        : x( x ), y( y ), w( w ), h( h ), c( c ), i( i ), j( j )
    {} // constructor
}; // struct

std::vector<rect> getCallOverview( std::shared_ptr<DBConSingle> pConnection, std::shared_ptr<Pack> pPack,
                                   int64_t iRunId, double dMinScore, int64_t iX, int64_t iY, uint64_t uiW, uint64_t uiH,
                                   uint64_t uiMaxW, uint64_t uiMaxH, uint32_t uiGiveUpFactor )
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
    size_t uiStartContigX = pPack->uiSequenceIdForPosition( uiX );
    size_t uiEndContigX = pPack->uiSequenceIdForPosition( uiX + uiW );
    size_t uiStartContigY = pPack->uiSequenceIdForPosition( uiY );
    size_t uiEndContigY = pPack->uiSequenceIdForPosition( uiY + uiH );
    std::vector<rect> vRet;
    for( size_t uiContigX = uiStartContigX; uiContigX <= uiEndContigX; uiContigX++ )
    {
        for( size_t uiContigY = uiStartContigY; uiContigY <= uiEndContigY; uiContigY++ )
        {
            auto uiStartX = std::max( uiX, (uint32_t)pPack->startOfSequenceWithId( uiContigX ) );
            auto uiEndX = std::min( uiX + uiW, (uint64_t)pPack->endOfSequenceWithId( uiContigX ) );
            auto uiStartY = std::max( uiY, (uint32_t)pPack->startOfSequenceWithId( uiContigY ) );
            auto uiEndY = std::min( uiY + uiH, (uint64_t)pPack->endOfSequenceWithId( uiContigY ) );
            uint32_t uiNumW = ( uint32_t )( uiEndX - uiStartX ) / (uint32_t)uiMaxW + (uint32_t)1;
            uint32_t uiNumH = ( uint32_t )( uiEndY - uiStartY ) / (uint32_t)uiMaxH + (uint32_t)1;
            double dW = ( (double)( uiEndX - uiStartX ) ) / (double)uiNumW;
            double dH = ( (double)( uiEndY - uiStartY ) ) / (double)uiNumH;
            if( dW * uiGiveUpFactor < uiW )
                continue;
            if( dH * uiGiveUpFactor < uiH )
                continue;
            for( size_t uiI = 0; uiI < uiNumW; uiI++ )
                for( size_t uiJ = 0; uiJ < uiNumH; uiJ++ )
                {
                    int64_t uiInnerX = ( int64_t )( uiI * dW + uiStartX );
                    int64_t uiInnerY = ( int64_t )( uiJ * dH + uiStartY );
                    auto c = getCallOverviewArea( pConnection, pPack, iRunId, dMinScore, uiInnerX, uiInnerY,
                                                  (uint32_t)dW + 1, (uint32_t)dH + 1 );
                    if( c > 0 )
                        vRet.emplace_back( (uint32_t)uiInnerX, (uint32_t)uiInnerY, (uint32_t)dW, (uint32_t)dH, c,
                                           (uint32_t)uiContigX, (uint32_t)uiContigY );
                } // for
        } // for
    } // for
    return vRet;
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
        .def( "add_score_index", &SvCallTable<DBConSingle>::addScoreIndex );

    py::class_<SvJumpRunTable<DBConSingle>, std::shared_ptr<SvJumpRunTable<DBConSingle>>>( rxPyModuleId,
                                                                                           "JumpRunTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) );

    py::class_<SvJumpTable<DBConSingle>, std::shared_ptr<SvJumpTable<DBConSingle>>>( rxPyModuleId, "SvJumpTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "create_indices", &SvJumpTable<DBConSingle>::createIndices )
        .def( "num_jumps", &SvJumpTable<DBConSingle>::numJumps );

    py::class_<ReadTable<DBConSingle>, std::shared_ptr<ReadTable<DBConSingle>>>( rxPyModuleId, "ReadTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "get_read", &ReadTable<DBConSingle>::getRead );

    py::class_<SvCallerRunTable<DBConSingle>, std::shared_ptr<SvCallerRunTable<DBConSingle>>>( rxPyModuleId,
                                                                                               "SvCallerRunTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "getIds", &SvCallerRunTable<DBConSingle>::getIds )
        .def( "getName", &SvCallerRunTable<DBConSingle>::getName )
        .def( "getDesc", &SvCallerRunTable<DBConSingle>::getDesc )
        .def( "jump_run_id", &SvCallerRunTable<DBConSingle>::getSvJumpRunId )
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
    rxPyModuleId.def( "get_call_overview", &getCallOverview );
    rxPyModuleId.def( "get_call_overview_area", &getCallOverviewArea );

    rxPyModuleId.def( "combine_overlapping_calls", &combineOverlappingCalls<DBConSingle> );

    exportSvCallInserter( rxPyModuleId );
    exportCallsFromDb( rxPyModuleId );
    exportSvJump( rxPyModuleId );
    exportSvJumpInserter( rxPyModuleId );
    exportNucSeqSql( rxPyModuleId );
    exportReadInserter( rxPyModuleId );
} // function
#endif // WITH_PYTHON
