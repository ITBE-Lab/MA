#include "container/sv_db/svDb.h"
#include "module/combineOverlappingCalls.h"
#include <pybind11/stl.h>

// include classes that implement sql queries
#include "container/sv_db/query_objects/callInserter.h"
#include "container/sv_db/query_objects/fetchCalls.h"
#include "container/sv_db/query_objects/fetchRuns.h"
#include "container/sv_db/query_objects/fetchSvJump.h"
#include "container/sv_db/query_objects/jumpInserter.h"
#include "container/sv_db/query_objects/nucSeqSql.h"
#include "container/sv_db/query_objects/readInserter.h"

using namespace libMA;

uint32_t getCallOverviewArea( std::shared_ptr<SV_DB> pDb, std::shared_ptr<Pack> pPack, int64_t iRunId, double dMinScore,
                              int64_t iX, int64_t iY, uint64_t uiW, uint64_t uiH )
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
    CppSQLiteExtQueryStatement<uint32_t> xQuery( *pDb->pDatabase,
                                                 "SELECT COUNT(*) "
                                                 "FROM sv_call_table, sv_call_r_tree "
                                                 "WHERE sv_call_table.id == sv_call_r_tree.id "
                                                 "AND sv_call_r_tree.run_id_b >= ? " // dim 1
                                                 "AND sv_call_r_tree.run_id_a <= ? " // dim 1
                                                 "AND sv_call_r_tree.maxX >= ? " // dim 2
                                                 "AND sv_call_r_tree.minX <= ? " // dim 2
                                                 "AND sv_call_r_tree.maxY >= ? " // dim 3
                                                 "AND sv_call_r_tree.minY <= ? " // dim 3
                                                 "AND (supporting_nt*1.0)/coverage >= ? " );
    return xQuery.scalar( iRunId, iRunId, uiX, uiX + (uint32_t)uiW, uiY, uiY + (uint32_t)uiH, dMinScore );
} // function

uint32_t getNumJumpsInArea( std::shared_ptr<SV_DB> pDb, std::shared_ptr<Pack> pPack, int64_t iRunId, int64_t iX,
                            int64_t iY, uint64_t uiW, uint64_t uiH )
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
    CppSQLiteExtQueryStatement<uint32_t> xQuery( *pDb->pDatabase,
                                                 "SELECT COUNT(*) "
                                                 "FROM sv_jump_table "
                                                 "WHERE sv_jump_run_id == ? "
                                                 "AND ( (from_pos >= ? AND from_pos <= ?) OR from_pos == ? ) "
                                                 "AND ( (to_pos >= ? AND to_pos <= ?) OR to_pos == ? ) "
                                                 "ORDER BY sort_pos_start" );
    return xQuery.scalar( iRunId, uiX, uiX + (uint32_t)uiW, std::numeric_limits<uint32_t>::max( ), uiY,
                          uiY + (uint32_t)uiH, std::numeric_limits<uint32_t>::max( ) );
} // function

struct rect
{
    uint32_t x, y, w, h, c, i, j;
    rect( uint32_t x, uint32_t y, uint32_t w, uint32_t h, uint32_t c, uint32_t i, uint32_t j )
        : x( x ), y( y ), w( w ), h( h ), c( c ), i( i ), j( j )
    {} // constructor
}; // struct

std::vector<rect> getCallOverview( std::shared_ptr<SV_DB> pDb, std::shared_ptr<Pack> pPack, int64_t iRunId,
                                   double dMinScore, int64_t iX, int64_t iY, uint64_t uiW, uint64_t uiH,
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
            uint32_t uiNumW = ( uiEndX - uiStartX ) / uiMaxW + 1;
            uint32_t uiNumH = ( uiEndY - uiStartY ) / uiMaxH + 1;
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
                    auto c = getCallOverviewArea( pDb, pPack, iRunId, dMinScore, uiInnerX, uiInnerY, (uint32_t)dW + 1,
                                                  (uint32_t)dH + 1 );
                    if( c > 0 )
                        vRet.emplace_back( uiInnerX, uiInnerY, (uint32_t)dW, (uint32_t)dH, c, uiContigX, uiContigY );
                } // for
        } // for
    } // for
    return vRet;
} // function

#ifdef WITH_PYTHON

void exportSoCDbWriter( py::module& rxPyModuleId )
{

    // export the SV_DB class
    py::class_<SV_DB, std::shared_ptr<SV_DB>>( rxPyModuleId, "SV_DB" ) //
        .def( py::init<std::string, std::string>( ) )
        .def( "delete_run", &SV_DB::deleteRun )
        .def( "get_run_id", &SV_DB::getRunId )
        .def( "get_call_area", &SV_DB::getCallArea )
        .def( "get_num_overlaps_between_calls", &SV_DB::getNumOverlapsBetweenCalls )
        .def( "get_blur_on_overlaps_between_calls", &SV_DB::getBlurOnOverlapsBetweenCalls )
        .def( "get_num_invalid_calls", &SV_DB::getNumInvalidCalls )
        .def( "add_score_index", &SV_DB::addScoreIndex )
        .def( "get_num_calls", &SV_DB::getNumCalls )
        .def( "get_run_name", &SV_DB::getRunName )
        .def( "get_run_desc", &SV_DB::getRunDesc )
        .def( "get_run_date", &SV_DB::getRunDate )
        .def( "get_run_jump_id", &SV_DB::getRunJumpId )
        .def( "get_num_runs", &SV_DB::getNumRuns )
        .def( "get_max_score", &SV_DB::getMaxScore )
        .def( "get_min_score", &SV_DB::getMinScore )
        .def( "get_num_nts", &SV_DB::getNumNts )
        .def( "run_exists", &SV_DB::runExists )
        .def( "name_exists", &SV_DB::nameExists )
        .def( "set_num_threads", &SV_DB::setNumThreads )
        .def( "create_jump_indices", &SV_DB::createJumpIndices )
        .def( "reconstruct_sequenced_genome", &SV_DB::reconstructSequencedGenome )
        .def( "newest_unique_runs", &SV_DB::getNewestUniqueRuns )
        .def( "update_coverage", &SV_DB::updateCoverage )
        .def( "insert_sv_caller_run", &SV_DB::insertSvCallerRun )
        .def( "insert_sv_jump_run", &SV_DB::insertSvJumpRun )
        .def( "get_read", &SV_DB::getRead )
        .def( "num_jumps", &SV_DB::numJumps );

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

    rxPyModuleId.def( "combine_overlapping_calls", &combineOverlappingCalls );

    exportSvCallInserter( rxPyModuleId );
    exportCallsFromDb( rxPyModuleId );
    exportRunsFromDb( rxPyModuleId );
    exportSvJump( rxPyModuleId );
    exportSvJumpInserter( rxPyModuleId );
    exportNucSeqSql( rxPyModuleId );
    exportReadInserter( rxPyModuleId );
} // function
#endif
