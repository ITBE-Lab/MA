#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#include "util/execution-context.h"
#include "util/export.h"
#include <cstdlib>
#include <iostream>

using namespace libMA;

std::shared_ptr<NucSeq> randomNucSeq( size_t uiLen )
{
    auto pRet = std::make_shared<NucSeq>( );
    pRet->vReserveMemory( uiLen );
    for( size_t i = 0; i < uiLen; i++ )
        pRet->push_back( (uint8_t)(std::rand( ) % 4) );
    return pRet;
} // function

int main( void )
{
    auto pPack = makePledge<Pack>( );
    pPack->get( )->vAppendSequence( "chr1", "chr1-desc", *randomNucSeq( 65536 ) );
    auto pFmIndex = makePledge<FMIndex>( pPack->get() );

    auto pQueryVec = std::make_shared<ContainerVector<std::shared_ptr<NucSeq>>>( );
    for( size_t i = 0; i < 1000; i++ )
        pQueryVec->push_back( randomNucSeq( 1000 ) );
    for( size_t i = 0; i < 100; i++ )
        pQueryVec->push_back( randomNucSeq( 10000 ) );

    ParameterSetManager xParameters;
    xParameters.getSelected()->xSearchInversions->set(true);

    auto vGraphSinks = setUpCompGraph(
        xParameters, pPack, pFmIndex, promiseMe( std::make_shared<StaticSplitter<NucSeq>>( xParameters, pQueryVec ) ),
        std::make_shared<Collector<NucSeq, ContainerVector<std::shared_ptr<Alignment>>, Pack>>( xParameters ),
        xParameters.pGeneralParameterSet->piNumberOfThreads->get( ) );

    BasePledge::simultaneousGet( vGraphSinks );

    return EXIT_SUCCESS;
} /// main function