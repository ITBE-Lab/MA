/**
 * @file export.cpp
 * @author Markus Schmidt
 */
#include "util/export.h"

using namespace libMA;

#ifdef WITH_PYTHON
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
    exportModule( );
    exportFM_index( );
    exportSequence( );
    exportBinarySeeding( );
    exportPack( );
    exportIntervalTree( );
    exportExceptions( );
    exportSeed( );
    exportAlignment( );
    exportLinesweep( );
    exportNeedlemanWunsch( );
    exportStripOfConsideration( );
    exportExtractAllSeeds( );
    exportExecOnVector( );
    exportFileReader( );
    exportFileWriter( );
    exportMappingQuality( );
    exportPairedReads( );
    exportSplitter( );
    exportSMW( );
    exportSW_GPU( );
    exportSoC( );
    exportOtherSeeding( );
    defaults::exportDefaults( );
} // function

#endif


std::vector<std::shared_ptr<BasePledge>> libMA::setUpCompGraph( std::shared_ptr<Pledge<Pack, false>> pPack,
                                                                std::shared_ptr<Pledge<FMIndex, false>>
                                                                    pFMDIndex,
                                                                std::shared_ptr<Pledge<NucSeq, true>>
                                                                    pQueries,
                                                                std::shared_ptr<WriterModule>
                                                                    pWriter,
                                                                unsigned int uiThreads )
{
    auto pLock = std::make_shared<Lock<NucSeq>>( );
    auto pSeeding = std::make_shared<BinarySeeding>( );
    auto pSOC = std::make_shared<StripOfConsideration>( );
    auto pHarmonization = std::make_shared<Harmonization>( );
    auto pDP = std::make_shared<NeedlemanWunsch>( );
    auto pMappingQual = std::make_shared<MappingQuality>( );
    std::vector<std::shared_ptr<BasePledge>> aRet;
    for( unsigned int i = 0; i < uiThreads; i++ )
    {
        auto pQuery = promiseMe( pLock, pQueries );
        auto pUnlock = std::make_shared<UnLock<Container>>( pQuery ); // might require the correct pledge...
        auto pSeeds = promiseMe( pSeeding, pFMDIndex, pQuery );
        auto pSOCs = promiseMe( pSOC, pSeeds, pQuery, pPack, pFMDIndex );
        auto pHarmonized = promiseMe( pHarmonization, pSOCs, pQuery );
        auto pAlignments = promiseMe( pDP, pHarmonized, pQuery, pPack );
        auto pAlignmentsWQuality = promiseMe( pMappingQual, pQuery, pAlignments );
        auto pEmptyContainer = promiseMe( pWriter, pQuery, pAlignmentsWQuality, pPack );
        auto pUnlockResult = promiseMe( pUnlock, pEmptyContainer );
        aRet.push_back( pUnlockResult );
    } // for
    return aRet;
} // function

#if 0
std::vector<std::shared_ptr<Pledge>> setUpCompGraph( std::shared_ptr<Pledge> pPack,
                                                     std::shared_ptr<Pledge>
                                                         pFMDIndex,
                                                     std::shared_ptr<Pledge>
                                                         pQueries,
                                                     std::vector<std::shared_ptr<Module>>& vOut,
                                                     unsigned int uiThreads )
{
    // setup all modules
    // modules required for any alignment
    std::shared_ptr<Module> pLockQuery( new Lock( std::shared_ptr<Container>( new NucSeq( ) ) ) );
    auto pSeeding = std::make_shared<BinarySeeding>( );
    auto pSOC = std::make_shared<StripOfConsideration>( );
    std::shared_ptr<LinearLineSweep> pCouple( new LinearLineSweep( ) );
    // we only want to report the best alignment
    std::shared_ptr<Module> pDoOptimal(
        new ExecOnVec( std::shared_ptr<Module>( new NeedlemanWunsch( ) ), true, 0 ) );
    std::shared_ptr<MappingQuality> pMapping( new MappingQuality( ) );


    auto pDummyQuery = std::make_shared<Pledge>( std::make_shared<NucSeq>( ) );
    pDummyQuery->set( std::make_shared<NucSeq>( ) );

    // setup the computational graph
    std::vector<std::shared_ptr<Pledge>> aRet;
    for( unsigned int i = 0; i < uiThreads; i++ )
    {
        // lock the query in for this subgraph
        std::shared_ptr<Pledge> pQuery =
            Module::promiseMe( pLockQuery, std::vector<std::shared_ptr<Pledge>>{pQueries} );
        // one more module to unlock the locked-in query
        // this requires the pledge from the lock
        // therefore we have to create one module for each subgraph
        std::shared_ptr<Module> pUnLock( new UnLock( pQuery ) );
        // the seeding stage
        std::shared_ptr<Pledge> pSeeds =
            Module::promiseMe( pSeeding, std::vector<std::shared_ptr<Pledge>>{pFMDIndex, pQuery} );
        // the filtering stage
        std::shared_ptr<Pledge> pSOCs = Module::promiseMe(
            pSOC, std::vector<std::shared_ptr<Pledge>>{pSeeds, pQuery, pPack, pFMDIndex} );
        // the coupling stage
        std::shared_ptr<Pledge> pCoupled =
            Module::promiseMe( pCouple, std::vector<std::shared_ptr<Pledge>>{pSOCs, pQuery} );
        if( defaults::bFindMode )
        {
            // write the output to a file
            assert( vOut.size( ) == 1 );
            std::shared_ptr<Pledge> pNil = Module::promiseMe(
                vOut[ 0 ], std::vector<std::shared_ptr<Pledge>>{pCoupled, pPack} );
            // unlock the query so that this subgraph can be executed multiple times
            std::shared_ptr<Pledge> pRet =
                Module::promiseMe( pUnLock, std::vector<std::shared_ptr<Pledge>>{pNil} );
            // save the
            aRet.push_back( pRet );
        } // if
        else
        {
            // the optimal matching stage
            std::shared_ptr<Pledge> pOptimal = Module::promiseMe(
                pDoOptimal, std::vector<std::shared_ptr<Pledge>>{pCoupled, pQuery, pPack} );
            // assign a mapping quality
            std::shared_ptr<Pledge> pAlignments = Module::promiseMe(
                pMapping, std::vector<std::shared_ptr<Pledge>>{pQuery, pOptimal} );
            // write the output to a file
            if( vOut.size( ) == 1 )
            {
                std::shared_ptr<Pledge> pNil = Module::promiseMe(
                    vOut[ 0 ],
                    std::vector<std::shared_ptr<Pledge>>{pQuery, pDummyQuery, pAlignments, pPack} );
                // unlock the query so that this subgraph can be executed multiple times
                std::shared_ptr<Pledge> pRet =
                    Module::promiseMe( pUnLock, std::vector<std::shared_ptr<Pledge>>{pNil} );
                // save the
                aRet.push_back( pRet );
            } // if
            else
            {
                assert( i < vOut.size( ) );
                std::shared_ptr<Pledge> pNil = Module::promiseMe(
                    vOut[ i ],
                    std::vector<std::shared_ptr<Pledge>>{pQuery, pDummyQuery, pAlignments, pPack} );
                // unlock the query so that this subgraph can be executed multiple times
                std::shared_ptr<Pledge> pRet =
                    Module::promiseMe( pUnLock, std::vector<std::shared_ptr<Pledge>>{pNil} );
                // save the pledge
                aRet.push_back( pRet );
            } // else
        } // else
    } // for

    return aRet;
} // function

std::vector<std::shared_ptr<Pledge>> setUpCompGraphPaired( std::shared_ptr<Pledge> pPack,
                                                           std::shared_ptr<Pledge>
                                                               pFMDIndex,
                                                           std::shared_ptr<Pledge>
                                                               pQueries,
                                                           std::vector<std::shared_ptr<libMA::Module>>& vOut,
                                                           unsigned int uiThreads )
{
    // setup all modules
    // modules required for any alignment
    std::shared_ptr<Module> pLockQueries( new Lock( std::make_shared<ContainerVector>( ) ) );
    std::shared_ptr<Module> pLock( new Lock( std::make_shared<NucSeq>( ) ) );
    auto pSplitter = std::make_shared<ReadSplitter>( );
    auto pSeeding = std::make_shared<BinarySeeding>( );
    auto pSOC = std::make_shared<StripOfConsideration>( );
    std::shared_ptr<LinearLineSweep> pCouple( new LinearLineSweep( ) );
    // we only want to report the best alignment
    std::shared_ptr<Module> pDoOptimal( new ExecOnVec( std::shared_ptr<Module>( new NeedlemanWunsch( ) ), true, 0 ) );
    std::shared_ptr<MappingQuality> pMapping( new MappingQuality( ) );
    // modules for the paired alignment
    std::shared_ptr<Module> pPaired( new PairedReads( ) );


    // setup the computational graph
    std::vector<std::shared_ptr<Pledge>> aRet;
    for( unsigned int i = 0; i < uiThreads; i++ )
    {
        // lock the query in for this subgraph
        std::shared_ptr<Pledge> pQueriesLocked =
            Module::promiseMe( pLockQueries, std::vector<std::shared_ptr<Pledge>>{pQueries} );
        std::shared_ptr<Pledge> pQuery =
            Module::promiseMe( pSplitter, std::vector<std::shared_ptr<Pledge>>{pQueriesLocked} );
        std::shared_ptr<Pledge> pQuery1 = Module::promiseMe( pLock, std::vector<std::shared_ptr<Pledge>>{pQuery} );
        std::shared_ptr<Pledge> pQuery2 = Module::promiseMe( pLock, std::vector<std::shared_ptr<Pledge>>{pQuery} );
        // one more module to unlock the locked-in query
        // this requires the pledge from the lock
        // therefore we have to create one module for each subgraph
        std::shared_ptr<Module> pUnLock( new UnLock( pQueriesLocked ) );
        std::shared_ptr<Module> pUnLock2( new UnLock( pQuery1 ) );
        std::shared_ptr<Module> pUnLock3( new UnLock( pQuery2 ) );
        // the seeding stage
        std::shared_ptr<Pledge> pSeeds1 =
            Module::promiseMe( pSeeding, std::vector<std::shared_ptr<Pledge>>{pFMDIndex, pQuery1} );
        std::shared_ptr<Pledge> pSeeds2 =
            Module::promiseMe( pSeeding, std::vector<std::shared_ptr<Pledge>>{pFMDIndex, pQuery2} );
        // the filtering stage
        std::shared_ptr<Pledge> pSOCs1 =
            Module::promiseMe( pSOC, std::vector<std::shared_ptr<Pledge>>{pSeeds1, pQuery1, pPack, pFMDIndex} );
        std::shared_ptr<Pledge> pSOCs2 =
            Module::promiseMe( pSOC, std::vector<std::shared_ptr<Pledge>>{pSeeds2, pQuery2, pPack, pFMDIndex} );
        // the coupling stage
        std::shared_ptr<Pledge> pCoupled1 =
            Module::promiseMe( pCouple, std::vector<std::shared_ptr<Pledge>>{pSOCs1, pQuery1} );
        std::shared_ptr<Pledge> pCoupled2 =
            Module::promiseMe( pCouple, std::vector<std::shared_ptr<Pledge>>{pSOCs2, pQuery2} );
        // the optimal matching stage
        std::shared_ptr<Pledge> pOptimal1 =
            Module::promiseMe( pDoOptimal, std::vector<std::shared_ptr<Pledge>>{pCoupled1, pQuery1, pPack} );
        std::shared_ptr<Pledge> pOptimal2 =
            Module::promiseMe( pDoOptimal, std::vector<std::shared_ptr<Pledge>>{pCoupled2, pQuery2, pPack} );
        // assign a mapping quality
        std::shared_ptr<Pledge> pAlignment1 =
            Module::promiseMe( pMapping, std::vector<std::shared_ptr<Pledge>>{pQuery1, pOptimal1} );
        std::shared_ptr<Pledge> pAlignment2 =
            Module::promiseMe( pMapping, std::vector<std::shared_ptr<Pledge>>{pQuery2, pOptimal2} );
        std::shared_ptr<Pledge> pAlignments =
            Module::promiseMe( pPaired, std::vector<std::shared_ptr<Pledge>>{pAlignment1, pAlignment2, pPack} );
        // write the output to a file
        if( vOut.size( ) == 1 )
        {
            // write the output to a file
            std::shared_ptr<Pledge> pNil = Module::promiseMe(
                vOut[ 0 ], std::vector<std::shared_ptr<Pledge>>{pQuery1, pQuery2, pAlignments, pPack} );
            // unlock the query so that this subgraph can be executed multiple times
            std::shared_ptr<Pledge> pRet = Module::promiseMe( pUnLock, std::vector<std::shared_ptr<Pledge>>{pNil} );
            // unlock the query so that this subgraph can be executed multiple times
            pRet = Module::promiseMe( pUnLock2, std::vector<std::shared_ptr<Pledge>>{pRet} );
            // unlock the query so that this subgraph can be executed multiple times
            pRet = Module::promiseMe( pUnLock3, std::vector<std::shared_ptr<Pledge>>{pRet} );
            // save the
            aRet.push_back( pRet );
        } // if
        else
        {
            assert( i < vOut.size( ) );
            // write the output to a file
            std::shared_ptr<Pledge> pNil = Module::promiseMe(
                vOut[ i ], std::vector<std::shared_ptr<Pledge>>{pQuery1, pQuery2, pAlignments, pPack} );
            // unlock the query so that this subgraph can be executed multiple times
            std::shared_ptr<Pledge> pRet = Module::promiseMe( pUnLock, std::vector<std::shared_ptr<Pledge>>{pNil} );
            // unlock the query so that this subgraph can be executed multiple times
            pRet = Module::promiseMe( pUnLock2, std::vector<std::shared_ptr<Pledge>>{pRet} );
            // unlock the query so that this subgraph can be executed multiple times
            pRet = Module::promiseMe( pUnLock3, std::vector<std::shared_ptr<Pledge>>{pRet} );
            // save the pledge
            aRet.push_back( pRet );
        } // else
    } // for

    return aRet;
} // function
#endif