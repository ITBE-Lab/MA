/**
 * @file export.cpp
 * @author Markus Schmidt
 */
#include "util/export.h"
#include "util/execution-context.h"
#include "util/parameter.h"

using namespace libMA;


#ifdef WITH_PYTHON

/*
 * The exposing of the ParameterSetManager and it's components has been reworked in SV branch
 * -> here only some quickfixes, so that python does not immedeately segfault if someone tries to use this...
 */

template <typename TP_VALUE> void exportAlignerParameter( py::module& rxPyModuleId, std::string sName )
{
    py::class_< //
        AlignerParameter<TP_VALUE>, //
        AlignerParameterBase, //
        std::shared_ptr<AlignerParameter<TP_VALUE>> //
        >( rxPyModuleId, sName.c_str( ) ) //
        .def( py::init<const std::string, const std::string, const size_t, const std::string, const TP_VALUE>() )
        .def( "set", &AlignerParameter<TP_VALUE>::set ) //
        .def( "get", &AlignerParameter<TP_VALUE>::get_py );
} // function

/**
 * Parameter set manager export is here, so that we do not need to include pybind11 in parameter.h
 */
void exportParameter( py::module& rxPyModuleId )
{
    py::class_<AlignerParameterBase, std::shared_ptr<AlignerParameterBase>>( rxPyModuleId, "AlignerParameterBase" ) //
        .def( "mirror", &AlignerParameterBase::mirror ) //
        .def_readonly( "name", &AlignerParameterBase::sName ) //
        .def_readonly( "description", &AlignerParameterBase::sDescription );

    exportAlignerParameter<int>( rxPyModuleId, "AlignerParameterInt" );
    exportAlignerParameter<short>( rxPyModuleId, "AlignerParameterShort" );
    exportAlignerParameter<bool>( rxPyModuleId, "AlignerParameterBool" );
    exportAlignerParameter<double>( rxPyModuleId, "AlignerParameterDouble" );
    exportAlignerParameter<float>( rxPyModuleId, "AlignerParameterFloat" );
    py::class_<AlignerParameter<AlignerParameterBase::ChoicesType>, AlignerParameterBase, std::shared_ptr<AlignerParameter<AlignerParameterBase::ChoicesType>>>(
        rxPyModuleId, "AlignerParameterChoice" ) //
        .def( "set", &AlignerParameter<AlignerParameterBase::ChoicesType>::set ) //
        .def( "get", &AlignerParameter<AlignerParameterBase::ChoicesType>::get_py );
    // exportAlignerParameter<fs::path>( rxPyModuleId, "AlignerParameterFilePath" ); @todo

    // Export ParameterSetBase Class
    py::class_<ParameterSetBase, std::shared_ptr<ParameterSetBase>>( rxPyModuleId, "ParameterSetBase" ) //
        .def( "mirror", &ParameterSetBase::mirror ) //
        .def( "register_parameter", &ParameterSetBase::registerParameter ) //
        .def( "unregister_parameter", &ParameterSetBase::unregisterParameter ) //
        .def( "by_name", &ParameterSetBase::byName ) //
        .def( "by_short", &ParameterSetBase::byShort );
    

    // Export Presetting Class
    py::class_<Presetting, ParameterSetBase, std::shared_ptr<Presetting>>( rxPyModuleId, "Presetting" ) //
        .def( py::init<std::string>( ) ) //
        .def_readonly( "name", &Presetting::sName );

    // Export Presetting Class
    py::class_<GeneralParameter, ParameterSetBase, std::shared_ptr<GeneralParameter>>( rxPyModuleId,
                                                                                       "GeneralSettings" ) //
        .def( py::init<>( ) );

    // Export ParameterSetManager Class
    py::class_<ParameterSetManager>( rxPyModuleId, "ParameterSetManager" ) //
        .def( py::init<>( ) ) //
        .def_readwrite( "global_settings", &ParameterSetManager::pGlobalParameterSet ) //
        .def( "get", &ParameterSetManager::get )
        .def( "add_setting", &ParameterSetManager::addSetting )
        .def( "get_num_threads", &ParameterSetManager::getNumThreads )
        .def( "set_selected", &ParameterSetManager::setSelected )
        .def( "get_selected", &ParameterSetManager::getSelected_py )
        .def( "by_name", &ParameterSetManager::byName )
        .def( "by_short", &ParameterSetManager::byShort );
} // function

void exportExecutionContext( py::module& rxPyModuleId )
{
    py::class_<GenomeManager>( rxPyModuleId, "GenomeManager" ) //
        .def( "load_genome", &GenomeManager::loadGenome );
    py::class_<ReadsManager>( rxPyModuleId, "ReadsManager" ) //
        .def_readwrite( "primary_queries", &ReadsManager::vsPrimaryQueryFullFileName ) //
        .def_readwrite( "mate_queries", &ReadsManager::vsMateQueryFullFileName );
    py::class_<OutputManager>( rxPyModuleId, "OutputManager" );

    py::class_<ExecutionContext>( rxPyModuleId, "ExecutionContext" ) //
        .def( py::init<>( ) ) //
        .def( "do_align", &ExecutionContext::doAlignCallbackLess ) //
        // def_readonly required since xParameterSetManager is non-copyable
        .def_readonly( "parameter_set_manager", &ExecutionContext::xParameterSetManager ) //
        .def_readwrite( "genome_manager", &ExecutionContext::xGenomeManager ) //
        .def_readwrite( "reads_manager", &ExecutionContext::xReadsManager ) //
        .def_readwrite( "output_manager", &ExecutionContext::xOutputManager );
} // function



/**
 * @brief pybind11 translator for AnnotatedException
 */
void translator( AnnotatedException const& x )
{
    PyErr_SetString( PyExc_RuntimeError, x.what( ) );
}
/**
 * @brief pybind11 exporter for AnnotatedException
 * @details
 * @todo once AnnotatedException is unused remove this export
 */
void exportExceptions( py::module& rxPyModuleId )
{
    py::register_exception<AnnotatedException>( rxPyModuleId, "AnnotatedException" );
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
    exportSVJump( libMaModule );
    exportSoCDbWriter( libMaModule );
    exportSoC( libMaModule );
    exportOtherSeeding( libMaModule );
    exportExecutionContext( libMaModule );
    exportHashMapSeeding( libMaModule );
    exportSvJumpsFromSeeds( libMaModule );
    exportSweepSvJump( libMaModule );
    exportConnectorPatternFilter( libMaModule );
    exportMinimizerIndex( libMaModule );
} // function

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
    auto pCast = std::make_shared<Cast<SuffixArrayInterface, FMIndex>>( rParameters );

    // create the graph
    std::vector<std::shared_ptr<BasePledge>> aRet;
    for( unsigned int i = 0; i < uiThreads; i++ )
    {
        auto pQuery = promiseMe( pLock, pQueries );
        auto pSeeds = promiseMe( pSeeding, promiseMe(pCast, pFMDIndex), pQuery );
        auto pSOCs = promiseMe( pSOC, pSeeds, pQuery, pPack, pFMDIndex );
        auto pHarmonized = promiseMe( pHarmonization, pSOCs, pQuery, pFMDIndex );
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
    auto pCast = std::make_shared<Cast<SuffixArrayInterface, FMIndex>>( rParameters );

    // create the graph
    std::vector<std::shared_ptr<BasePledge>> aRet;
    for( unsigned int i = 0; i < uiThreads; i++ )
    {
        auto pQueryTuple = promiseMe( pLock, pQueries );
        auto pQueryA = promiseMe( pGetFirst, pQueryTuple );
        auto pQueryB = promiseMe( pGetSecond, pQueryTuple );
        auto pSeedsA = promiseMe( pSeeding, promiseMe(pCast, pFMDIndex), pQueryA );
        auto pSeedsB = promiseMe( pSeeding, promiseMe(pCast, pFMDIndex), pQueryB );
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
