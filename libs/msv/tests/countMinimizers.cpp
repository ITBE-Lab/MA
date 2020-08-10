#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#include "msv/module/count_k_mers.h"
#include <iostream>
#include <map>
#include <time.h>

using namespace libMA;
using namespace libMSV;


#define MAX_VAL (2ull << 56)
int main( void )
{
    srand( static_cast<unsigned int>( time( NULL ) ) );


    std::map<uint64_t, size_t> vTestSet;

    auto pCnt = std::make_shared<HashCounter>( );

    for( size_t uiI = 0; uiI < 100; uiI++ )
    {
        uint64_t uiHash = std::rand( ) % MAX_VAL;
        size_t uiCnt = ( std::rand( ) % 3 ) + 1;
        vTestSet[ uiHash ] = uiCnt;
        std::cout << "inserting: " << uiHash << " " << uiCnt << std::endl;
    }

    for( auto& xPair : vTestSet )
        for( size_t uiI = 0; uiI < xPair.second; uiI++ )
            pCnt->addHash( xPair.first );

    pCnt->iterate(
        [&]( uint64_t uiHash, size_t uiCnt ) { std::cout << "have: " << uiHash << " " << uiCnt << std::endl; } );

    // check that all inserted sequences actually exist
    for( auto& xPair : vTestSet )
        if( pCnt->isUnique( xPair.first ) )
            assert( false );

    // check that no new sequences are in the HashCounter
    for( size_t uiI = 0; uiI < 100; uiI++ )
    {
        uint64_t uiHash = std::rand( ) % MAX_VAL;
        if( vTestSet.count( uiHash ) == 0 )
            assert( pCnt->isUnique( uiHash ) );
    } // for

    // check that the backcomputing of NucSeqs works and that all sequences are extracted correctly
    pCnt->iterate( [&]( uint64_t uiHash, size_t uiCnt ) {
        auto xIt = vTestSet.find( uiHash );
        assert( xIt->second == uiCnt );
        vTestSet.erase( xIt );
    } );
    assert( vTestSet.size( ) == 0 );

    return EXIT_SUCCESS;
} /// main function
