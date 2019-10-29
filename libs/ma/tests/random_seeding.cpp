#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#include "util/execution-context.h"
#include "util/export.h"
#include <cstdlib>
#include <iostream>

#if _MSC_VER

// @todo the template class Join does nut work with msvc
int main( void )
{
    return EXIT_SUCCESS;
} /// main function

#elif __GNUC__

using namespace libMA;

std::shared_ptr<NucSeq> randomNucSeq( size_t uiLen )
{
    auto pRet = std::make_shared<NucSeq>( );
    pRet->vReserveMemory( uiLen );
    for( size_t i = 0; i < uiLen; i++ )
        pRet->push_back( ( uint8_t )( std::rand( ) % 4 ) );
    return pRet;
} // function

class SegVecChecker : public Module<Container, false, SegmentVector, NucSeq, Pack, FMIndex>
{
  public:
    SegVecChecker( const ParameterSetManager& rParameters )
    {} // constructor

    virtual std::shared_ptr<Container> execute( std::shared_ptr<SegmentVector> pIn, std::shared_ptr<NucSeq> pQuery,
                                                std::shared_ptr<Pack> pRef, std::shared_ptr<FMIndex> pFmIndex )
    {
        pIn->extractSeeds( pFmIndex, 1000, 16 )->confirmSeedPositions( pQuery, pRef, true );
        return std::make_shared<Container>( );
    } // method
}; // class

#ifdef WITH_ZLIB
class SeedsChecker : public Module<Container, false, Seeds, NucSeq, Pack, FMIndex>
{
  public:
    SeedsChecker( const ParameterSetManager& rParameters )
    {} // constructor

    virtual std::shared_ptr<Container> execute( std::shared_ptr<Seeds> pIn, std::shared_ptr<NucSeq> pQuery,
                                                std::shared_ptr<Pack> pRef, std::shared_ptr<FMIndex> pFmIndex )
    {
        pIn->confirmSeedPositions( pQuery, pRef, false );
        return std::make_shared<Container>( );
    } // method
}; // class
#endif

template class Join<Container, Container
#ifdef WITH_ZLIB
                    ,
                    Container
#endif
                    >;

int main( void )
{
    auto pPack = makePledge<Pack>( );
    pPack->get( )->vAppendSequence( "chr1", "chr1-desc", *randomNucSeq( 65536 ) );
    auto pFmIndex = makePledge<FMIndex>( pPack->get( ) );

    auto pQueryVec = std::make_shared<ContainerVector<std::shared_ptr<NucSeq>>>( );
    for( size_t i = 0; i < 1000; i++ )
        pQueryVec->push_back( randomNucSeq( 100 ) );
    for( size_t i = 0; i < 100; i++ )
        pQueryVec->push_back( randomNucSeq( 1000 ) );

    ParameterSetManager xParameters;

    // set up the modules
    auto pLock = std::make_shared<Lock<NucSeq>>( xParameters );
    auto pFMDIndexSeeding1 = std::make_shared<BinarySeeding>( xParameters );
    xParameters.getSelected( )->xSeedingTechnique->set( 1 ); // use SMEM seeding
    auto pFMDIndexSeeding2 = std::make_shared<BinarySeeding>( xParameters );
    auto pChecker = std::make_shared<SegVecChecker>( xParameters );
#ifdef WITH_ZLIB
    auto pMMIndex =
        makePledge<minimizer::Index>( xParameters, pPack->get( )->contigSeqs( ), pPack->get( )->contigNames( ) );
    auto pMMISeeding = std::make_shared<MinimizerSeeding>( xParameters );
    auto pChecker2 = std::make_shared<SeedsChecker>( xParameters );
#endif
    auto pJoin = std::make_shared<Join<Container, Container, Container>>( xParameters );

    auto pQueries = promiseMe( std::make_shared<StaticSplitter<NucSeq>>( xParameters, pQueryVec ) );

    // create the graph
    std::vector<std::shared_ptr<BasePledge>> vGraphSinks;
    for( int i = 0; i < xParameters.xGlobalParameterSet.piNumberOfThreads->get( ); i++ )
    {
        auto pQuery = promiseMe( pLock, pQueries );

        // check the max spanning seeds
        auto pSeeds = promiseMe( pFMDIndexSeeding1, pFmIndex, pQuery );
        auto pEmpty = promiseMe( pChecker, pSeeds, pQuery, pPack, pFmIndex );

        // check the SMEMs
        auto pSeeds2 = promiseMe( pFMDIndexSeeding2, pFmIndex, pQuery );
        auto pEmpty2 = promiseMe( pChecker, pSeeds2, pQuery, pPack, pFmIndex );

#ifdef WITH_ZLIB
        // check the minimizers
        auto pSeeds3 = promiseMe( pMMISeeding, pMMIndex, pQuery, pPack );
        auto pEmpty3 = promiseMe( pChecker2, pSeeds3, pQuery, pPack, pFmIndex );
#endif

        // wait for both tasks to be completed
        auto pJoined = promiseMe( pJoin, pEmpty, pEmpty2
#ifdef WITH_ZLIB
                                  ,
                                  pEmpty3
#endif
        );
        // unlock the query once completed
        auto pUnlockResult = promiseMe( std::make_shared<UnLock<Container>>( xParameters, pQuery ), pJoined );

        // add this task to todo list
        vGraphSinks.push_back( pUnlockResult );
    } // for

    BasePledge::simultaneousGet( vGraphSinks );

    return EXIT_SUCCESS;
} /// main function

#endif