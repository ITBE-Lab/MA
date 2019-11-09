/**
 * @file export.cpp
 * @author Markus Schmidt
 */
#include "util/export.h"
#include "util/parameter.h"

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
    exportSegment( );
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
    exportSoC( );
    exportOtherSeeding( );
    defaults::exportDefaults( );
} // function

#else


/*
 * The exposing of the ParameterSetManager and it's components has been reworked in SV branch
 * -> here only some quickfixes, so that python does not immedeately segfault if someone tries to use this...
 */

template <typename TP_VALUE> void exportAlignerParameter( py::module& rxPyModuleId, std::string sName )
{
    py::class_<AlignerParameter<TP_VALUE>, AlignerParameterBase, std::shared_ptr<AlignerParameter<TP_VALUE>>>(
        rxPyModuleId, sName.c_str( ) ) //
        .def( "set", &AlignerParameter<TP_VALUE>::set ) //
        .def( "get", &AlignerParameter<TP_VALUE>::get_py );
} // function

/**
 * Parameter set manager export is here, so that we do not need to include pybind11 in parameter.h
 */
void exportParameter( py::module& rxPyModuleId )
{
    py::class_<AlignerParameterBase, std::shared_ptr<AlignerParameterBase>>( rxPyModuleId, "AlignerParameterBase" ) //
        .def_readonly( "name", &AlignerParameterBase::sName ) //
        .def_readonly( "description", &AlignerParameterBase::sDescription );

    exportAlignerParameter<int>( rxPyModuleId, "AlignerParameterInt" );
    exportAlignerParameter<short>( rxPyModuleId, "AlignerParameterShort" );
    exportAlignerParameter<bool>( rxPyModuleId, "AlignerParameterBool" );
    exportAlignerParameter<float>( rxPyModuleId, "AlignerParameterFloat" );
    // exportAlignerParameter<?>(rxPyModuleId, "AlignerParameterFilePath" );
    // exportAlignerParameter<?>(rxPyModuleId, "AlignerParameterChoice" );
    
    py::class_<AlignerParameter<AlignerParameterBase::ChoicesType>, AlignerParameterBase, std::shared_ptr<AlignerParameter<AlignerParameterBase::ChoicesType>>>(
        rxPyModuleId, "AlignerParameterChoice" ) //
        .def( "set", &AlignerParameter<AlignerParameterBase::ChoicesType>::set ) //
        .def( "get", &AlignerParameter<AlignerParameterBase::ChoicesType>::get_py );

    // Export Presetting Class
    py::class_<Presetting>( rxPyModuleId, "Presetting" ) //
        .def( py::init<>( ) ); //
    //.def( "__setitem__", &Presetting::byName )
    //.def( "__getitem__", &Presetting::byName );

    // Export ParameterSetManager Class
    py::class_<ParameterSetManager>( rxPyModuleId, "ParameterSetManager" ) //
        .def( py::init<>( ) ) //
        .def( "get", &ParameterSetManager::get )
        .def( "by_name", &ParameterSetManager::byName )
        .def( "by_short", &ParameterSetManager::byShort )
        .def( "set_selected", &ParameterSetManager::setSelected );
    //.def( "get_selected", &ParameterSetManager::getSelected_py )
} // function

PYBIND11_MODULE( libMA, libMaModule )
{
    DEBUG_3( std::cout.setf( std::ios::unitbuf ); )
    exportParameter( libMaModule );
    exportContainer( libMaModule );
    exportModuleClass( libMaModule );
    exportFM_index( libMaModule );
    exportNucSeq( libMaModule );
    exportBinarySeeding( libMaModule );
    exportPack( libMaModule );
    exportSegment( libMaModule );
    exportExceptions( libMaModule );
    exportSeed( libMaModule );
    exportAlignment( libMaModule );
    exportHarmonization( libMaModule );
    exportNeedlemanWunsch( libMaModule );
    exportStripOfConsideration( libMaModule );
    exportFileReader( libMaModule );
    exportFileWriter( libMaModule );
    exportMappingQuality( libMaModule );
    exportPairedReads( libMaModule );
    exportSplitter( libMaModule );
    exportSmallInversions( libMaModule );
#ifdef WITH_POSTGRES
    exportDBWriter( libMaModule );
#endif
    exportSoC( libMaModule );
    exportOtherSeeding( libMaModule );
    defaults::exportDefaults( libMaModule );
#ifdef WITH_ZLIB
    exportMinimizerIndex( libMaModule );
    exportMinimizerSeeding( libMaModule );
#endif
} // function

#endif

#endif


std::vector<std::shared_ptr<BasePledge>> libMA::setUpCompGraph( const ParameterSetManager& rParameters,
                                                                std::shared_ptr<Pledge<Pack>>
                                                                    pPack,
                                                                std::shared_ptr<Pledge<FMIndex>>
                                                                    pFMDIndex,
                                                                std::shared_ptr<Pledge<NucSeq, true>>
                                                                    pQueries,
                                                                std::shared_ptr<TP_WRITER>
                                                                    pWriter,
                                                                unsigned int uiThreads )
{
    // set up the modules
    auto pLock = std::make_shared<Lock<NucSeq>>( rParameters );
    auto pSeeding = std::make_shared<BinarySeeding>( rParameters );
    auto pSOC = std::make_shared<StripOfConsideration>( rParameters );
    auto pHarmonization = std::make_shared<Harmonization>( rParameters );
    auto pDP = std::make_shared<NeedlemanWunsch>( rParameters );
    auto pMappingQual = std::make_shared<MappingQuality>( rParameters );
    auto pSmallInversions = std::make_shared<SmallInversions>( rParameters );

    // create the graph
    std::vector<std::shared_ptr<BasePledge>> aRet;
    for( unsigned int i = 0; i < uiThreads; i++ )
    {
        auto pQuery = promiseMe( pLock, pQueries );
        auto pSeeds = promiseMe( pSeeding, pFMDIndex, pQuery );
        auto pSOCs = promiseMe( pSOC, pSeeds, pQuery, pPack, pFMDIndex );
        auto pHarmonized = promiseMe( pHarmonization, pSOCs, pQuery );
        auto pAlignments = promiseMe( pDP, pHarmonized, pQuery, pPack );
        auto pAlignmentsWQuality = promiseMe( pMappingQual, pQuery, pAlignments );
        if( rParameters.getSelected( )->xSearchInversions->get( ) )
        {
            auto pAlignmentsWInv = promiseMe( pSmallInversions, pAlignmentsWQuality, pQuery, pPack );
            auto pEmptyContainer = promiseMe( pWriter, pQuery, pAlignmentsWInv, pPack );
            auto pUnlockResult =
                promiseMe( std::make_shared<UnLock<Container>>( rParameters, pQuery ), pEmptyContainer );
            aRet.push_back( pUnlockResult );
        } // if
        else
        {
            auto pEmptyContainer = promiseMe( pWriter, pQuery, pAlignmentsWQuality, pPack );
            auto pUnlockResult =
                promiseMe( std::make_shared<UnLock<Container>>( rParameters, pQuery ), pEmptyContainer );
            aRet.push_back( pUnlockResult );
        } // else
    } // for
    return aRet;
} // function

#ifdef WITH_ZLIB
std::vector<std::shared_ptr<BasePledge>> libMA::setUpCompGraph( const ParameterSetManager& rParameters,
                                                                std::shared_ptr<Pledge<Pack>>
                                                                    pPack,
                                                                std::shared_ptr<Pledge<minimizer::Index>>
                                                                    pIndex,
                                                                std::shared_ptr<Pledge<NucSeq, true>>
                                                                    pQueries,
                                                                std::shared_ptr<TP_WRITER>
                                                                    pWriter,
                                                                unsigned int uiThreads )
{
    // set up the modules
    auto pLock = std::make_shared<Lock<NucSeq>>( rParameters );
    auto pSeeding = std::make_shared<MinimizerSeeding>( rParameters );
    auto pSOC = std::make_shared<MinimizerStripOfConsideration>( rParameters );
    auto pHarmonization = std::make_shared<Harmonization>( rParameters );
    auto pDP = std::make_shared<NeedlemanWunsch>( rParameters );
    auto pMappingQual = std::make_shared<MappingQuality>( rParameters );
    auto pSmallInversions = std::make_shared<SmallInversions>( rParameters );

    // create the graph
    std::vector<std::shared_ptr<BasePledge>> aRet;
    for( unsigned int i = 0; i < uiThreads; i++ )
    {
        auto pQuery = promiseMe( pLock, pQueries );
        auto pSeeds = promiseMe( pSeeding, pIndex, pQuery, pPack );
        auto pSOCs = promiseMe( pSOC, pSeeds, pQuery, pPack );
        auto pHarmonized = promiseMe( pHarmonization, pSOCs, pQuery );
        auto pAlignments = promiseMe( pDP, pHarmonized, pQuery, pPack );
        auto pAlignmentsWQuality = promiseMe( pMappingQual, pQuery, pAlignments );
        if( rParameters.getSelected( )->xSearchInversions->get( ) )
        {
            auto pAlignmentsWInv = promiseMe( pSmallInversions, pAlignmentsWQuality, pQuery, pPack );
            auto pEmptyContainer = promiseMe( pWriter, pQuery, pAlignmentsWInv, pPack );
            auto pUnlockResult =
                promiseMe( std::make_shared<UnLock<Container>>( rParameters, pQuery ), pEmptyContainer );
            aRet.push_back( pUnlockResult );
        } // if
        else
        {
            auto pEmptyContainer = promiseMe( pWriter, pQuery, pAlignmentsWQuality, pPack );
            auto pUnlockResult =
                promiseMe( std::make_shared<UnLock<Container>>( rParameters, pQuery ), pEmptyContainer );
            aRet.push_back( pUnlockResult );
        } // else
    } // for
    return aRet;
} // function
#endif

std::vector<std::shared_ptr<BasePledge>> libMA::setUpCompGraphPaired( const ParameterSetManager& rParameters,
                                                                      std::shared_ptr<Pledge<Pack>>
                                                                          pPack,
                                                                      std::shared_ptr<Pledge<FMIndex>>
                                                                          pFMDIndex,
                                                                      std::shared_ptr<Pledge<TP_PAIRED_READS, true>>
                                                                          pQueries,
                                                                      std::shared_ptr<TP_PAIRED_WRITER>
                                                                          pWriter,
                                                                      unsigned int uiThreads )
{
    // set up the modules
    auto pLock = std::make_shared<Lock<TP_PAIRED_READS>>( rParameters );
    auto pGetFirst = std::make_shared<TupleGet<TP_PAIRED_READS, 0>>( rParameters );
    auto pGetSecond = std::make_shared<TupleGet<TP_PAIRED_READS, 1>>( rParameters );
    auto pSeeding = std::make_shared<BinarySeeding>( rParameters );
    auto pSOC = std::make_shared<StripOfConsideration>( rParameters );
    auto pHarmonization = std::make_shared<Harmonization>( rParameters );
    auto pDP = std::make_shared<NeedlemanWunsch>( rParameters );
    auto pMappingQual = std::make_shared<MappingQuality>( rParameters );
    auto pSmallInversions = std::make_shared<SmallInversions>( rParameters );
    auto pPairedReads = std::make_shared<PairedReads>( rParameters );

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
        auto pHarmonizedA = promiseMe( pHarmonization, pSOCsA, pQueryA );
        auto pHarmonizedB = promiseMe( pHarmonization, pSOCsB, pQueryB );
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
            auto pUnlockResult =
                promiseMe( std::make_shared<UnLock<Container>>( rParameters, pQueryTuple ), pEmptyContainer );
            aRet.push_back( pUnlockResult );
        } // if
        else
        {
            auto pAlignmentsWQuality =
                promiseMe( pPairedReads, pQueryA, pQueryB, pAlignmentsWQualityA, pAlignmentsWQualityB, pPack );
            auto pEmptyContainer = promiseMe( pWriter, pQueryA, pQueryB, pAlignmentsWQuality, pPack );
            auto pUnlockResult =
                promiseMe( std::make_shared<UnLock<Container>>( rParameters, pQueryTuple ), pEmptyContainer );
            aRet.push_back( pUnlockResult );
        } // else
    } // for
    return aRet;
} // function
