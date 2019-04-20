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
    for(size_t i = 0; i < 100; i++)
    {
        auto pPack = std::make_shared<Pack>( );
        for(size_t i = 0; i <= ((size_t)std::rand( ) % 5); i++)
            pPack->vAppendSequence( "chr1", "chr1-desc", *randomNucSeq( (uint8_t)(std::rand( ) % 30000) + 1 ) );
        auto pFmIndex = std::make_shared<FMIndex>( pPack );

        if(!pFmIndex->test(pPack, 1000))
            return EXIT_FAILURE;
    } // for

    return EXIT_SUCCESS;
} /// main function