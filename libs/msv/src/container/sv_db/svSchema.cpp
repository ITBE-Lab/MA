#include "msv/container/sv_db/svSchema.h"
#include "ms/container/container.h"
#include "ms/container/sv_db/py_db_conf.h"
#include "msv/module/combineOverlappingCalls.h"

// include classes that implement sql queries
#include "ms/util/pybind11.h"
#include "msv/container/sv_db/query_objects/callInserter.h"
#include "msv/container/sv_db/query_objects/fetchCalls.h"
#include "msv/container/sv_db/query_objects/fetchSvJump.h"
#include "msv/container/sv_db/query_objects/jumpInserter.h"
#include "msv/container/sv_db/query_objects/kMerInserter.h"
#include "msv/container/sv_db/query_objects/nucSeqSql.h"
#include "msv/container/sv_db/query_objects/readInserter.h"
#include "pybind11/stl.h"

using namespace libMSV;


template <typename DBCon>
uint32_t getNumJumpsInArea( std::shared_ptr<DBCon> pConnection, std::shared_ptr<Pack> pPack, int64_t iRunId, int64_t iX,
                            int64_t iY, uint64_t uiW, uint64_t uiH, uint64_t uiLimit )
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


    auto xRectangle = WKBUint64Rectangle( geom::Rectangle<nucSeqIndex>( uiX, uiY, uiW, uiH ) );

    SQLQuery<DBConSingle, uint64_t> xQuery( pConnection,
                                            "SELECT COUNT(*) FROM ( "
                                            // requires subquery here so that limit actually optimizes the count
                                            "   SELECT id "
                                            "   FROM sv_jump_table "
                                            "   WHERE sv_jump_run_id = ? "
                                            "   AND " ST_INTERSCTS "(rectangle, ST_PolyFromWKB(?, 0)) "
                                            "   LIMIT ? "
                                            ") AS tmp_table " );
    // FIXME: Don't use scalar anymore!
    return xQuery.scalar( iRunId, xRectangle, uiLimit );
} // function


#ifdef WITH_PYTHON

void exportSoCDbWriter( libMS::SubmoduleOrganizer& xOrganizer )
{
    py::class_<SvCallTable<DBConSingle>, std::shared_ptr<SvCallTable<DBConSingle>>>( xOrganizer.util( ), "SvCallTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "reconstruct_sequenced_genome", &SvCallTable<DBConSingle>::reconstructSequencedGenome )
        .def( "calls_to_seeds", &SvCallTable<DBConSingle>::callsToSeeds )
        .def( "calls_to_seeds_by_id", &SvCallTable<DBConSingle>::callsToSeedsById )
        .def( "calls_to_seeds_by_id_auto", &SvCallTable<DBConSingle>::callsToSeedsByIdAutoStart )
        .def( "reconstruct_sequenced_genome_from_seeds",
              &SvCallTable<DBConSingle>::reconstructSequencedGenomeFromSeeds )
        .def( "num_calls", &SvCallTable<DBConSingle>::numCalls_py )
        .def( "max_score", &SvCallTable<DBConSingle>::maxScore )
        .def( "min_score", &SvCallTable<DBConSingle>::minScore )
        .def( "call_area", &SvCallTable<DBConSingle>::callArea )
        .def( "drop_indices", &SvCallTable<DBConSingle>::dropIndices )
        .def( "filter_calls_with_high_score", &SvCallTable<DBConSingle>::filterCallsWithHighScore )
        .def( "gen_indices", &SvCallTable<DBConSingle>::genIndices )
        .def( "insert_call", &SvCallTable<DBConSingle>::insertCall );

    using X = SvCallTableAnalyzer<DBCon, false>;
    py::class_<X, std::shared_ptr<X>>( xOrganizer.util( ), "SvCallTableAnalyzer" )
        .def( py::init<std::shared_ptr<PoolContainer<DBCon>>>( ) )
        .def( "num_overlaps", &X::numOverlaps )
        .def( "num_invalid_calls", &X::numInvalidCalls )
        .def( "blur_on_overlaps", &X::blurOnOverlaps );

    py::class_<OneSidedCallsTable<DBConSingle>, std::shared_ptr<OneSidedCallsTable<DBConSingle>>>(
        xOrganizer.util( ), "OneSidedCallsTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "get_name", &OneSidedCallsTable<DBConSingle>::getMate )
        .def( "insert_calls", &OneSidedCallsTable<DBConSingle>::insertCalls );

    py::class_<KMerFilterTable<DBConSingle>, std::shared_ptr<KMerFilterTable<DBConSingle>>>( xOrganizer.util( ),
                                                                                             "KMerFilterTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "get_counter", &KMerFilterTable<DBConSingle>::getCounter )
        .def( "insert_counter_set", &KMerFilterTable<DBConSingle>::insert_counter_set );

    py::class_<HashFilterTable<DBConSingle>, std::shared_ptr<HashFilterTable<DBConSingle>>>( xOrganizer.util( ),
                                                                                             "HashFilterTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "get_counter", &HashFilterTable<DBConSingle>::getCounter )
        .def( "insert_counter_set", &HashFilterTable<DBConSingle>::insert_counter_set );

    py::class_<SvJumpRunTable<DBConSingle>, std::shared_ptr<SvJumpRunTable<DBConSingle>>>( xOrganizer.util( ),
                                                                                           "JumpRunTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) );
    py::class_<SequencerTable<DBConSingle>, std::shared_ptr<SequencerTable<DBConSingle>>>( xOrganizer.util( ),
                                                                                           "SequencerTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) );
    py::class_<CallDescTable<DBConSingle>, std::shared_ptr<CallDescTable<DBConSingle>>>( xOrganizer.util( ),
                                                                                         "CallDescTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "insert", &CallDescTable<DBConSingle>::insert_py )
        .def( "gen_index", &CallDescTable<DBConSingle>::genIndex )
        .def( "get_desc", &CallDescTable<DBConSingle>::getDesc );

    py::class_<SvJumpTable<DBConSingle>, std::shared_ptr<SvJumpTable<DBConSingle>>>( xOrganizer.util( ), "SvJumpTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "create_indices", &SvJumpTable<DBConSingle>::createIndices )
        .def( "drop_indices", &SvJumpTable<DBConSingle>::dropIndices )
        .def( "num_jumps", &SvJumpTable<DBConSingle>::numJumps );

    py::class_<ReadTable<DBConSingle>, std::shared_ptr<ReadTable<DBConSingle>>>( xOrganizer.util( ), "ReadTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "get_read", &ReadTable<DBConSingle>::getRead )
        .def( "get_read_id", &ReadTable<DBConSingle>::getReadId )
        .def( "get_seq_id", &ReadTable<DBConSingle>::getSeqId )
        .def( "get_used_reads", &ReadTable<DBConSingle>::getUsedReads );

    py::class_<ReadRangeTable<DBConSingle>, std::shared_ptr<ReadRangeTable<DBConSingle>>>( xOrganizer.util( ),
                                                                                           "ReadRangeTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( py::init<std::shared_ptr<DBConSingle>, bool>( ) )
        .def( "insert", &ReadRangeTable<DBConSingle>::insertAlignment )
        .def( "insert", &ReadRangeTable<DBConSingle>::insertAlignmentId )
        .def( "coverage", &ReadRangeTable<DBConSingle>::coverage )
        .def( "drop_indices", &ReadRangeTable<DBConSingle>::dropIndices )
        .def( "gen_indices", &ReadRangeTable<DBConSingle>::genIndices );

    py::class_<ReadSelectionTable<DBConSingle>, std::shared_ptr<ReadSelectionTable<DBConSingle>>>( xOrganizer.util( ),
                                                                                           "ReadSelectionTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "insert_by_name", &ReadSelectionTable<DBConSingle>::insertReadByName );

    py::class_<SvCallSupportTable<DBConSingle>, std::shared_ptr<SvCallSupportTable<DBConSingle>>>(
        xOrganizer.util( ), "SvCallSupportTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) );

    py::class_<SvCallerRunTable<DBConSingle>, std::shared_ptr<SvCallerRunTable<DBConSingle>>>( xOrganizer.util( ),
                                                                                               "SvCallerRunTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "getIds", &SvCallerRunTable<DBConSingle>::getIds )
        .def( "getName", &SvCallerRunTable<DBConSingle>::getName )
        .def( "exists", &SvCallerRunTable<DBConSingle>::exists )
        .def( "getDesc", &SvCallerRunTable<DBConSingle>::getDesc )
        .def( "jump_run_id", &SvCallerRunTable<DBConSingle>::getSvJumpRunId )
        .def( "newest_unique_runs", &SvCallerRunTable<DBConSingle>::getNewestUnique )
        .def( "getDate", &SvCallerRunTable<DBConSingle>::getDate );

    py::class_<rect>( xOrganizer.util( ), "rect" ) //
        .def_readonly( "x", &rect::x )
        .def_readonly( "y", &rect::y )
        .def_readonly( "w", &rect::w )
        .def_readonly( "h", &rect::h )
        .def_readonly( "i", &rect::i )
        .def_readonly( "j", &rect::j )
        .def_readonly( "c", &rect::c );
    py::bind_vector<std::vector<rect>>( xOrganizer.util( ), "rectVector", "docstr" );
    xOrganizer.util( ).def( "get_num_jumps_in_area", &getNumJumpsInArea<DBConSingle> );
    xOrganizer.util( ).def( "get_call_overview", &getCallOverview<DBCon> );
    xOrganizer.util( ).def( "get_call_overview_area", &getCallOverviewArea<DBConSingle> );

    xOrganizer.util( ).def( "combine_overlapping_calls", &combineOverlappingCalls<DBCon> );

    py::class_<SQLDBInformer<DBConSingle>, std::shared_ptr<SQLDBInformer<DBConSingle>>>( xOrganizer.util( ),
                                                                                         "SQLDBInformer" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "get_all_schemas", &SQLDBInformer<DBConSingle>::getAllSchemas );

    exportSvCallInserter( xOrganizer );
    exportCallsFromDb( xOrganizer );
    exportSvJump( xOrganizer );
    exportSvJumpInserter( xOrganizer );
    exportNucSeqSql( xOrganizer );
    exportReadInserter( xOrganizer );
    exportKMerInserter( xOrganizer );
} // function
#endif // WITH_PYTHON
