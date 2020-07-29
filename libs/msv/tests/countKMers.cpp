#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#include "ma/container/nucSeq.h"
#include "msv/module/count_k_mers.h"
#include <iostream>
#include <map>
#include <time.h>

using namespace libMA;
using namespace libMSV;

std::shared_ptr<NucSeq> randomNucSeq( size_t uiLen, // length sequence
                                      size_t uiNMod = 0, // N start, statistically; 0 = no N at all
                                      size_t uiMaxSizeNseq = 150 ) // max size of N sequence

{
    auto pRet = std::make_shared<NucSeq>( );
    pRet->vReserveMemory( uiLen );
    int iNCounter = 0;
    for( size_t i = 0; i < uiLen; i++ )
    {
        if( uiNMod )
            if( iNCounter == 0 && std::rand( ) % uiNMod == 0 )
            {
                iNCounter = std::rand( ) % uiMaxSizeNseq;
            }
        if( iNCounter > 0 )
        {
            iNCounter--;
            pRet->push_back( 4 ); // put N's intentionally!
        }
        else
            pRet->push_back( ( uint8_t )( std::rand( ) % 4 ) );
    }
    return pRet;
} // function

int main( void )
{
    // srand( static_cast<unsigned int>( time( NULL ) ) );


    std::map<std::string, size_t> vTestSet;

    auto pCnt = std::make_shared<KMerCounter>( 1 );

    for( size_t uiI = 0; uiI < 1; uiI++ )
    {
        auto pNuc = randomNucSeq( 18 );
        size_t uiCnt = ( std::rand( ) % 3 ) + 1;
        vTestSet[ pNuc->toString( ) ] = uiCnt;
        std::cout << "inserting: " << pNuc->toString( ) << " " << uiCnt << std::endl;
        pNuc->vReverseAll( );
        pNuc->vSwitchAllBasePairsToComplement( );
        vTestSet[ pNuc->toString( ) ] = uiCnt;
        std::cout << "inserting: " << pNuc->toString( ) << " " << uiCnt << std::endl;
    }
    pCnt->iterate( [&]( const NucSeq& xSeq, size_t uiCnt ) {
        std::cout << "have: " << xSeq.toString( ) << " " << uiCnt << std::endl;
    } );
    for( auto& xPair : vTestSet )
        for( size_t uiI = 0; uiI < xPair.second; uiI++ )
            pCnt->addSequence( std::make_shared<NucSeq>( xPair.first ) );

    for( auto& xPair : vTestSet )
        assert( !pCnt->isUnique( std::make_shared<NucSeq>( xPair.first ) ) );

    for( size_t uiI = 0; uiI < 100; uiI++ )
    {
        auto pSeq = randomNucSeq( 18 );
        if( vTestSet.count( pSeq->toString( ) ) == 0 )
            assert( pCnt->isUnique( pSeq ) );
    } // for

    pCnt->iterate( [&]( const NucSeq& xSeq, size_t uiCnt ) {
        auto xIt = vTestSet.find( xSeq.toString( ) );
        assert( xIt->second == uiCnt * 2 );
        vTestSet.erase( xIt );
    } );
    assert( vTestSet.size( ) == 0 );

    return EXIT_SUCCESS;
} /// main function
