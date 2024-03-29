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
#include "msv/container/sv_db/query_objects/coverageUpdater.h"
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

    if( uiW == 0 && uiH == 0 )
        return 0;

    WKBUint64Rectangle xRectangle = geom::Rectangle<nucSeqIndex>( uiX, uiY, uiW, uiH );

    SQLQuery<DBConSingle, uint64_t> xQuery( pConnection,
                                            "SELECT COUNT(*) FROM ( "
                                            // requires subquery here so that limit actually optimizes the count
                                            "   SELECT id "
                                            "   FROM sv_jump_table "
                                            "   WHERE sv_jump_run_id = ? "
                                            "   AND " ST_INTERSCTS "(rectangle, ST_GeomFromWKB(?, 0)) "
                                            "   LIMIT ? "
                                            ") AS tmp_table " );

    auto xRet = xQuery.scalar( iRunId, xRectangle, uiLimit );
    return xRet;
} // function

std::shared_ptr<Pack> reconstructSequencedGenomeFromSeedsBorderd(
    std::vector<std::tuple<std::string, std::shared_ptr<Seeds>, std::vector<std::shared_ptr<NucSeq>>>>&
        vReconstructedSeeds,
    std::shared_ptr<Pack>
        pRef,
    nucSeqIndex uiContigBorder )
{
    auto pRet = std::make_shared<Pack>( );
    for( auto& xTup : vReconstructedSeeds )
    {
        NucSeq xCurrChrom;
        for( size_t uiI = 0; uiI < std::get<1>( xTup )->size( ); uiI++ )
        {
            auto xSeed = ( *std::get<1>( xTup ) )[ uiI ];

            if( xSeed.bOnForwStrand )
                pRef->vExtractSubsectionN( xSeed.start_ref( ), xSeed.end_ref( ), xCurrChrom, true );
            else
                pRef->vExtractSubsectionN( pRef->uiPositionToReverseStrand( xSeed.start_ref( ) ),
                                           pRef->uiPositionToReverseStrand( xSeed.start_ref( ) - xSeed.size( ) ),
                                           xCurrChrom,
                                           true );

            if( !std::get<2>( xTup ).empty( ) )
            {
                auto pNucSeq = std::get<2>( xTup )[ uiI ];
                if( pNucSeq != nullptr && pNucSeq->length( ) > 0 )
                    xCurrChrom.vAppend( pNucSeq->pxSequenceRef, pNucSeq->length( ) );
            } // if
        } // for

        pRet->vAppendSequence( std::get<0>( xTup ), "no_description", xCurrChrom );
    } // for
    return pRet;
} // method

std::shared_ptr<Pack> reconstructSequencedGenomeFromSeeds(
    std::vector<std::tuple<std::string, std::shared_ptr<Seeds>, std::vector<std::shared_ptr<NucSeq>>>>
        vReconstructedSeeds,
    std::shared_ptr<Pack>
        pRef )
{
    return reconstructSequencedGenomeFromSeedsBorderd( vReconstructedSeeds, pRef, 0 );
} // method


#ifdef WITH_PYTHON

void exportSoCDbWriter( libMS::SubmoduleOrganizer& xOrganizer )
{
    // @todo this should be a module...
    xOrganizer.util().def("reconstruct_sequenced_genome", reconstructSequencedGenomeFromSeeds);

    py::class_<SvCallTable<DBConSingle>, std::shared_ptr<SvCallTable<DBConSingle>>>( xOrganizer.util( ), "SvCallTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "calls_to_seeds_by_id", &SvCallTable<DBConSingle>::callsToSeedsById )
        .def( "num_calls", &SvCallTable<DBConSingle>::numCalls_py )
        .def( "copy_path", &SvCallTable<DBConSingle>::copyPath )
        .def( "max_score", &SvCallTable<DBConSingle>::maxScore )
        .def( "min_score", &SvCallTable<DBConSingle>::minScore )
        .def( "call_area", &SvCallTable<DBConSingle>::callArea )
        .def( "drop_indices", &SvCallTable<DBConSingle>::dropIndices )
        .def( "filter_calls_with_high_score", &SvCallTable<DBConSingle>::filterCallsWithHighScore )
        .def( "gen_indices", &SvCallTable<DBConSingle>::genIndices )
        .def( "extract_small_calls", &SvCallTable<DBConSingle>::extractSmallCalls )
        .def( "extract_large_calls", &SvCallTable<DBConSingle>::extractLargeCalls )
        .def( "insert_call", &SvCallTable<DBConSingle>::insertCall );

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


    py::class_<CoverageCollector, libMS::Container, std::shared_ptr<CoverageCollector>>( xOrganizer.util( ), "CoverageCollector" )
        .def( py::init<libMA::nucSeqIndex>( ) );

    py::class_<CoverageTable<DBConSingle>, std::shared_ptr<CoverageTable<DBConSingle>>>( xOrganizer.util( ),
                                                                                           "CoverageTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "get_coverage", &CoverageTable<DBConSingle>::getCoverage )
        .def( "drop_indices", &CoverageTable<DBConSingle>::dropIndices )
        .def( "create_indices", &CoverageTable<DBConSingle>::createIndices )
        .def( "init_coverage", &CoverageTable<DBConSingle>::initCoverage )
        .def( "init_coverage_from_col", &CoverageTable<DBConSingle>::initCoverageFromColl )
        ;

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
        .def( "read_name", &ReadTable<DBConSingle>::readName )
        .def( "get_used_reads", &ReadTable<DBConSingle>::getUsedReads );

    py::class_<ReadRangeTable<DBConSingle>, std::shared_ptr<ReadRangeTable<DBConSingle>>>( xOrganizer.util( ),
                                                                                           "ReadRangeTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( py::init<std::shared_ptr<DBConSingle>, bool>( ) )
        .def( "insert", &ReadRangeTable<DBConSingle>::insertAlignment )
        .def( "insert", &ReadRangeTable<DBConSingle>::insertAlignmentId )
        .def( "insert_range", &ReadRangeTable<DBConSingle>::insertRange )
        .def( "coverage", &ReadRangeTable<DBConSingle>::coverage )
        .def( "coverage", &ReadRangeTable<DBConSingle>::coveragePrim )
        .def( "drop_indices", &ReadRangeTable<DBConSingle>::dropIndices )
        .def( "gen_indices", &ReadRangeTable<DBConSingle>::genIndices );

    py::class_<ReadSelectionTable<DBConSingle>, std::shared_ptr<ReadSelectionTable<DBConSingle>>>(
        xOrganizer.util( ), "ReadSelectionTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "insert_by_name", &ReadSelectionTable<DBConSingle>::insertReadByName );

    py::class_<SvCallSupportTable<DBConSingle>, std::shared_ptr<SvCallSupportTable<DBConSingle>>>(
        xOrganizer.util( ), "SvCallSupportTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) );

    py::class_<SvCallerRunTable<DBConSingle>, std::shared_ptr<SvCallerRunTable<DBConSingle>>>( xOrganizer.util( ),
                                                                                               "SvCallerRunTable" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "getIds", &SvCallerRunTable<DBConSingle>::getIds )
        .def( "getId", &SvCallerRunTable<DBConSingle>::getId )
        .def( "hasName", &SvCallerRunTable<DBConSingle>::hasName )
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
    xOrganizer.util( ).def( "merge_dummy_calls", &mergeDummyCalls<DBCon> );

    py::class_<SQLDBInformer<DBConSingle>, std::shared_ptr<SQLDBInformer<DBConSingle>>>( xOrganizer.util( ),
                                                                                         "SQLDBInformer" )
        .def( py::init<std::shared_ptr<DBConSingle>>( ) )
        .def( "get_all_schemas", &SQLDBInformer<DBConSingle>::getAllSchemas );

    exportSvCallInserter( xOrganizer );
    exportCallsFromDb( xOrganizer );
    exportSvJump( xOrganizer );
    exportSvJumpInserter( xOrganizer );
    exportCoverageUpdater( xOrganizer );
    exportNucSeqSql( xOrganizer );
    exportReadInserter( xOrganizer );
    exportKMerInserter( xOrganizer );
} // function
#endif // WITH_PYTHON
