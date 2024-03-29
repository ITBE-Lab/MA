/**
 * @file export.cpp
 * @author Markus Schmidt
 */
#include "ma/util/export.h"
#include "ma/container/minimizer_index.h"
#include "ma/container/soc.h"
#include "ma/util/execution-context.h"
#include "ms/module/splitter.h"
#include "ms/util/export.h"
#include "ms/util/parameter.h"

using namespace libMA;
using namespace libMS;


#ifdef WITH_PYTHON

void exportExecutionContext( libMS::SubmoduleOrganizer& xOrganizer )
{
    py::class_<GenomeManager>( xOrganizer.util( ), "GenomeManager" ) //
        .def( "load_genome", &GenomeManager::loadGenome );
    py::class_<ReadsManager>( xOrganizer.util( ), "ReadsManager" ) //
        .def_readwrite( "primary_queries", &ReadsManager::vsPrimaryQueryFullFileName ) //
        .def_readwrite( "mate_queries", &ReadsManager::vsMateQueryFullFileName );
    py::class_<OutputManager>( xOrganizer.util( ), "OutputManager" );

    py::class_<ExecutionContext>( xOrganizer.util( ), "ExecutionContext" ) //
        .def( py::init<>( ) ) //
        .def( "do_align", &ExecutionContext::doAlignCallbackLess ) //
        // def_readonly required since xParameterSetManager is non-copyable
        .def_readonly( "parameter_set_manager", &ExecutionContext::xParameterSetManager ) //
        .def_readwrite( "genome_manager", &ExecutionContext::xGenomeManager ) //
        .def_readwrite( "reads_manager", &ExecutionContext::xReadsManager ) //
        .def_readwrite( "output_manager", &ExecutionContext::xOutputManager );
} // function

PYBIND11_MODULE( libMA, libMaModule )
{
    py::module::import( "libMS" );
    libMS::SubmoduleOrganizer xOrganizer( libMaModule );

    exportFM_index( xOrganizer );
    exportNucSeq( xOrganizer );
    exportBinarySeeding( xOrganizer );
    exportPack( xOrganizer );
    exportSegment( xOrganizer );
    exportSeed( xOrganizer );
    exportAlignment( xOrganizer );
    exportHarmonization( xOrganizer );
    exportNeedlemanWunsch( xOrganizer );
    exportStripOfConsideration( xOrganizer );
    exportFileReader( xOrganizer );
    exportFileWriter( xOrganizer );
    exportMappingQuality( xOrganizer );
    exportPairedReads( xOrganizer );
    exportSplitter( xOrganizer );
    exportSmallInversions( xOrganizer );
    exportSoC( xOrganizer );
    exportOtherSeeding( xOrganizer );
    exportHashMapSeeding( xOrganizer );
#ifdef WITH_ZLIB
    exportMinimizerIndex( xOrganizer );
#endif
    exportExecutionContext( xOrganizer );
    exportSamFileReader( xOrganizer );
    exportCompareAlignments( xOrganizer );
    
#ifdef WITH_ZLIB
    exportMinimizerSeeding( xOrganizer );
#endif
} // function

#endif


std::vector<std::shared_ptr<BasePledge>> libMA::setUpCompGraph( const ParameterSetManager& rParameters,
                                                                std::shared_ptr<Pledge<Pack>>
                                                                    pPack,
                                                                std::shared_ptr<Pledge<FMIndex>>
                                                                    pFMDIndex,
                                                                std::shared_ptr<Pledge<FileStreamQueue, false>>
                                                                    pQueue,
                                                                std::shared_ptr<TP_WRITER>
                                                                    pWriter,
                                                                unsigned int uiThreads )
{
    // set up the modules
    auto pFileStreamPicker = std::make_shared<QueuePicker<FileStream>>( rParameters );
    auto pLock = std::make_shared<Lock<FileStream>>( rParameters );
    auto pFileReader = std::make_shared<FileReader>( rParameters );
    auto pFileStreamPlacer = std::make_shared<QueuePlacer<NucSeq, FileStream>>( rParameters );
    auto pSeeding = std::make_shared<BinarySeeding>( rParameters );
    auto pSOC = std::make_shared<StripOfConsideration>( rParameters );
    auto pHarmonization = std::make_shared<Harmonization>( rParameters );
    auto pDP = std::make_shared<NeedlemanWunsch>( rParameters );
    auto pMappingQual = std::make_shared<MappingQuality>( rParameters );
    auto pSmallInversions = std::make_shared<SmallInversions>( rParameters );
    auto pCast = std::make_shared<Cast<SuffixArrayInterface, FMIndex>>( rParameters );
    auto pProgressPrinter = std::make_shared<ProgressPrinter<FileStreamQueue>>( rParameters );

    // create the graph
    std::vector<std::shared_ptr<BasePledge>> aRet;
    BasePledge::parallelGraph( uiThreads, [&]( ) {
        auto pPickedFile = promiseMe( pFileStreamPicker, pQueue );
        auto pLockedFile = promiseMe( pLock, pPickedFile );
        auto pQuery_ = promiseMe( pFileReader, pLockedFile );
        auto pQuery = promiseMe( pFileStreamPlacer, pQuery_, pLockedFile, pQueue );
        auto pSeeds = promiseMe( pSeeding, promiseMe( pCast, pFMDIndex ), pQuery );
        auto pSOCs = promiseMe( pSOC, pSeeds, pQuery, pPack, pFMDIndex );
        auto pHarmonized = promiseMe( pHarmonization, pSOCs, pQuery, pFMDIndex );
        auto pAlignments = promiseMe( pDP, pHarmonized, pQuery, pPack );
        auto pAlignmentsWQuality = promiseMe( pMappingQual, pQuery, pAlignments );
        if( rParameters.getSelected( )->xSearchInversions->get( ) )
        {
            auto pAlignmentsWInv = promiseMe( pSmallInversions, pAlignmentsWQuality, pQuery, pPack );
            auto pEmptyContainer = promiseMe( pWriter, pQuery, pAlignmentsWInv, pPack );
            auto pEmptyContainer_ = promiseMe( pProgressPrinter, pEmptyContainer, pQueue );
            auto pUnlockResult =
                promiseMe( std::make_shared<UnLock<libMS::Container>>( rParameters, pQuery ), pEmptyContainer_ );
            aRet.push_back( pUnlockResult );
        } // if
        else
        {
            auto pEmptyContainer = promiseMe( pWriter, pQuery, pAlignmentsWQuality, pPack );
            auto pEmptyContainer_ = promiseMe( pProgressPrinter, pEmptyContainer, pQueue );
            auto pUnlockResult =
                promiseMe( std::make_shared<UnLock<libMS::Container>>( rParameters, pLockedFile ), pEmptyContainer_ );
            aRet.push_back( pUnlockResult );
        } // else
    } ); // parallelGraph
    return aRet;
} // function

std::vector<std::shared_ptr<BasePledge>>
libMA::setUpCompGraphPaired( const ParameterSetManager& rParameters,
                             std::shared_ptr<Pledge<Pack>>
                                 pPack,
                             std::shared_ptr<Pledge<FMIndex>>
                                 pFMDIndex,
                             std::shared_ptr<Pledge<PairedFileStreamQueue, false>>
                                 pQueue,
                             std::shared_ptr<TP_PAIRED_WRITER>
                                 pWriter,
                             unsigned int uiThreads )
{
    // set up the modules
    auto pFileStreamPicker = std::make_shared<QueuePicker<PairedFileStream>>( rParameters );
    auto pLock = std::make_shared<Lock<PairedFileStream>>( rParameters );
    auto pFileReader = std::make_shared<PairedFileReader>( rParameters );
    auto pFileStreamPlacer = std::make_shared<QueuePlacer<PairedReadsContainer, PairedFileStream>>( rParameters );
    auto pGetFirst = std::make_shared<TupleGet<PairedReadsContainer, 0>>( rParameters );
    auto pGetSecond = std::make_shared<TupleGet<PairedReadsContainer, 1>>( rParameters );
    auto pSeeding = std::make_shared<BinarySeeding>( rParameters );
    auto pSOC = std::make_shared<StripOfConsideration>( rParameters );
    auto pHarmonization = std::make_shared<Harmonization>( rParameters );
    auto pDP = std::make_shared<NeedlemanWunsch>( rParameters );
    auto pMappingQual = std::make_shared<MappingQuality>( rParameters );
    auto pSmallInversions = std::make_shared<SmallInversions>( rParameters );
    auto pPairedReads = std::make_shared<PairedReads>( rParameters );
    auto pCast = std::make_shared<Cast<SuffixArrayInterface, FMIndex>>( rParameters );
    auto pProgressPrinter = std::make_shared<ProgressPrinter<PairedFileStreamQueue>>( rParameters );

    // create the graph
    std::vector<std::shared_ptr<BasePledge>> aRet;
    BasePledge::parallelGraph( uiThreads, [&]( ) {
        auto pPickedFile = promiseMe( pFileStreamPicker, pQueue );
        auto pLockedFile = promiseMe( pLock, pPickedFile );
        auto pQueryTuple_ = promiseMe( pFileReader, pLockedFile );
        auto pQueryTuple = promiseMe( pFileStreamPlacer, pQueryTuple_, pLockedFile, pQueue );
        auto pQueryA = promiseMe( pGetFirst, pQueryTuple );
        auto pQueryB = promiseMe( pGetSecond, pQueryTuple );
        auto pSeedsA = promiseMe( pSeeding, promiseMe( pCast, pFMDIndex ), pQueryA );
        auto pSeedsB = promiseMe( pSeeding, promiseMe( pCast, pFMDIndex ), pQueryB );
        auto pSOCsA = promiseMe( pSOC, pSeedsA, pQueryA, pPack, pFMDIndex );
        auto pSOCsB = promiseMe( pSOC, pSeedsB, pQueryB, pPack, pFMDIndex );
        auto pHarmonizedA = promiseMe( pHarmonization, pSOCsA, pQueryA, pFMDIndex );
        auto pHarmonizedB = promiseMe( pHarmonization, pSOCsB, pQueryB, pFMDIndex );
        auto pAlignmentsA = promiseMe( pDP, pHarmonizedA, pQueryA, pPack );
        auto pAlignmentsB = promiseMe( pDP, pHarmonizedB, pQueryB, pPack );
        auto pAlignmentsWQualityA = promiseMe( pMappingQual, pQueryA, pAlignmentsA );
        auto pAlignmentsWQualityB = promiseMe( pMappingQual, pQueryB, pAlignmentsB );
        if( rParameters.getSelected( )->xSearchInversions->get( ) )
        {
            auto pAlignmentsWInvA = promiseMe( pSmallInversions, pAlignmentsWQualityA, pQueryA, pPack );
            auto pAlignmentsWInvB = promiseMe( pSmallInversions, pAlignmentsWQualityB, pQueryB, pPack );
            auto pAlignmentsWQuality =
                promiseMe( pPairedReads, pQueryA, pQueryB, pAlignmentsWInvA, pAlignmentsWInvB, pPack );
            auto pEmptyContainer = promiseMe( pWriter, pQueryA, pQueryB, pAlignmentsWQuality, pPack );
            auto pEmptyContainer_ = promiseMe( pProgressPrinter, pEmptyContainer, pQueue );
            auto pUnlockResult =
                promiseMe( std::make_shared<UnLock<libMS::Container>>( rParameters, pLockedFile ), pEmptyContainer_ );
            aRet.push_back( pUnlockResult );
        } // if
        else
        {
            auto pAlignmentsWQuality =
                promiseMe( pPairedReads, pQueryA, pQueryB, pAlignmentsWQualityA, pAlignmentsWQualityB, pPack );
            auto pEmptyContainer = promiseMe( pWriter, pQueryA, pQueryB, pAlignmentsWQuality, pPack );
            auto pEmptyContainer_ = promiseMe( pProgressPrinter, pEmptyContainer, pQueue );
            auto pUnlockResult =
                promiseMe( std::make_shared<UnLock<libMS::Container>>( rParameters, pLockedFile ), pEmptyContainer_ );
            aRet.push_back( pUnlockResult );
        } // else
    } ); // parallelGraph
    return aRet;
} // function
