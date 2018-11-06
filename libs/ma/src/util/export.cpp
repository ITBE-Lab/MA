/**
 * @file export.cpp
 * @author Markus Schmidt
 */
#include "util/export.h"

using namespace libMA;

#ifdef WITH_PYTHON

#ifdef BOOST_PYTHON
/**
 * @brief The boost-python main method.
 *
 * A main function that exposes all Containers and Modules to python
 * by calling the respective export*() functions.
 * The main method is generated using boost-python.
 */
BOOST_PYTHON_MODULE( libMA )
{
    DEBUG_3( std::cout.setf( std::ios::unitbuf ); )
    exportContainer( );
    exportModuleClass( );
    exportFM_index( );
    exportNucSeq( );
    exportBinarySeeding( );
    exportPack( );
    exportIntervalTree( );
    exportExceptions( );
    exportSeed( );
    exportAlignment( );
    exportHarmonization( );
    exportNeedlemanWunsch( );
    exportStripOfConsideration( );
    exportFileReader( );
    exportFileWriter( );
    exportMappingQuality( );
    exportPairedReads( );
    exportSplitter( );
#ifdef WITH_POSTGRES
    exportDBWriter( );
#endif
    exportSoCDbWriter();
    exportSoC( );
    exportOtherSeeding( );
    defaults::exportDefaults( );
} // function

#else

PYBIND11_MODULE(libMA, libMaModule) {
    DEBUG_3( std::cout.setf( std::ios::unitbuf ); )
    exportContainer( libMaModule );
    exportModuleClass( libMaModule );
    //exportFM_index( libMaModule );
    exportNucSeq( libMaModule );
    //exportBinarySeeding( libMaModule );
    //exportPack( libMaModule );
    //exportIntervalTree( libMaModule );
    //exportExceptions( libMaModule );
    //exportSeed( libMaModule );
    //exportAlignment( libMaModule );
    //exportHarmonization( libMaModule );
    //exportNeedlemanWunsch( libMaModule );
    //exportStripOfConsideration( libMaModule );
    //exportFileReader( libMaModule );
    //exportFileWriter( libMaModule );
    //exportMappingQuality( libMaModule );
    //exportPairedReads( libMaModule );
    //exportSplitter( libMaModule );
#ifdef WITH_POSTGRES
    //exportDBWriter( libMaModule );
#endif
    //exportSoCDbWriter( libMaModule );
    //exportSoC( libMaModule );
    //exportOtherSeeding( libMaModule );
    //defaults::exportDefaults( libMaModule );
} // function

#endif

#endif


std::vector<std::shared_ptr<BasePledge>> libMA::setUpCompGraph( std::shared_ptr<Pledge<Pack>> pPack,
                                                                std::shared_ptr<Pledge<FMIndex>>
                                                                    pFMDIndex,
                                                                std::shared_ptr<Pledge<NucSeq, true>>
                                                                    pQueries,
                                                                std::shared_ptr<TP_WRITER>
                                                                    pWriter,
                                                                unsigned int uiThreads )
{
    // set up the modules
    auto pLock = std::make_shared<Lock<NucSeq>>( );
    auto pSeeding = std::make_shared<BinarySeeding>( );
    auto pSOC = std::make_shared<StripOfConsideration>( );
    auto pHarmonization = std::make_shared<Harmonization>( );
    auto pDP = std::make_shared<NeedlemanWunsch>( );
    auto pMappingQual = std::make_shared<MappingQuality>( );

    // create the graph
    std::vector<std::shared_ptr<BasePledge>> aRet;
    for( unsigned int i = 0; i < uiThreads; i++ )
    {
        auto pQuery = promiseMe( pLock, pQueries );
        auto pSeeds = promiseMe( pSeeding, pFMDIndex, pQuery );
        auto pSOCs = promiseMe( pSOC, pSeeds, pQuery, pPack, pFMDIndex );
        auto pHarmonized = promiseMe( pHarmonization, pSOCs, pQuery, pFMDIndex );
        auto pAlignments = promiseMe( pDP, pHarmonized, pQuery, pPack );
        auto pAlignmentsWQuality = promiseMe( pMappingQual, pQuery, pAlignments );
        auto pEmptyContainer = promiseMe( pWriter, pQuery, pAlignmentsWQuality, pPack );
        auto pUnlockResult = promiseMe( std::make_shared<UnLock<Container>>( pQuery ), pEmptyContainer );
        aRet.push_back( pUnlockResult );
    } // for
    return aRet;
} // function

std::vector<std::shared_ptr<BasePledge>> libMA::setUpCompGraphPaired( std::shared_ptr<Pledge<Pack>> pPack,
                                                                      std::shared_ptr<Pledge<FMIndex>>
                                                                          pFMDIndex,
                                                                      std::shared_ptr<Pledge<TP_PAIRED_READS, true>>
                                                                          pQueries,
                                                                      std::shared_ptr<TP_PAIRED_WRITER>
                                                                          pWriter,
                                                                      unsigned int uiThreads )
{
    // set up the modules
    auto pLock = std::make_shared<Lock<TP_PAIRED_READS>>( );
    auto pGetFirst = std::make_shared<TupleGet<TP_PAIRED_READS, 0>>( );
    auto pGetSecond = std::make_shared<TupleGet<TP_PAIRED_READS, 1>>( );
    auto pSeeding = std::make_shared<BinarySeeding>( );
    auto pSOC = std::make_shared<StripOfConsideration>( );
    auto pHarmonization = std::make_shared<Harmonization>( );
    auto pDP = std::make_shared<NeedlemanWunsch>( );
    auto pMappingQual = std::make_shared<MappingQuality>( );
    auto pPairedReads = std::make_shared<PairedReads>( );

    // create the graph
    std::vector<std::shared_ptr<BasePledge>> aRet;
    for( unsigned int i = 0; i < uiThreads; i++ )
    {
        auto pQueryTuple = promiseMe( pLock, pQueries );
        auto pQueryA = promiseMe( pGetFirst, pQueryTuple );
        auto pQueryB = promiseMe( pGetSecond, pQueryTuple );
        auto pSeedsA = promiseMe( pSeeding, pFMDIndex, pQueryA );
        auto pSeedsB = promiseMe( pSeeding, pFMDIndex, pQueryB );
        auto pSOCsA = promiseMe( pSOC, pSeedsA, pQueryA, pPack, pFMDIndex );
        auto pSOCsB = promiseMe( pSOC, pSeedsB, pQueryB, pPack, pFMDIndex );
        auto pHarmonizedA = promiseMe( pHarmonization, pSOCsA, pQueryA, pFMDIndex );
        auto pHarmonizedB = promiseMe( pHarmonization, pSOCsB, pQueryB, pFMDIndex );
        auto pAlignmentsA = promiseMe( pDP, pHarmonizedA, pQueryA, pPack );
        auto pAlignmentsB = promiseMe( pDP, pHarmonizedB, pQueryB, pPack );
        auto pAlignmentsWQualityA = promiseMe( pMappingQual, pQueryA, pAlignmentsA );
        auto pAlignmentsWQualityB = promiseMe( pMappingQual, pQueryB, pAlignmentsB );
        auto pAlignmentsWQuality = promiseMe( pPairedReads, pAlignmentsA, pAlignmentsB, pPack );
        auto pEmptyContainer = promiseMe( pWriter, pQueryA, pQueryB, pAlignmentsWQuality, pPack );
        auto pUnlockResult = promiseMe( std::make_shared<UnLock<Container>>( pQueryTuple ), pEmptyContainer );
        aRet.push_back( pUnlockResult );
    } // for
    return aRet;
} // function
