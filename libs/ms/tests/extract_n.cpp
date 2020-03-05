#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#include "util/execution-context.h"
#include "util/export.h"
//#include <cstdlib>
#include <iostream>

using namespace libMA;

std::shared_ptr<NucSeq> randomNucSeq( size_t uiLen )
{
    auto pRet = std::make_shared<NucSeq>( );
    pRet->vReserveMemory( uiLen );
    for( size_t i = 0; i < uiLen; i++ )
        pRet->push_back( ( uint8_t )( std::rand( ) % 5 ) ); // put N's intentionally!
    return pRet;
} // function

int main( void )
{
    for( size_t i = 0; i < 100; i++ )
    {
        auto pPack = std::make_shared<Pack>( );
        auto pNucSeq = randomNucSeq( ( uint8_t )( std::rand( ) % 3000 ) + 20 );
        pPack->vAppendSequence( "chr1", "chr1-desc", *pNucSeq );

        NucSeq xSeq2;

        for( size_t i = 0; i < 100; i++ )
        {
            nucSeqIndex iStart = ( nucSeqIndex )( std::rand( ) % ( pNucSeq->length( ) - 10 ) + 10 );
            nucSeqIndex iEnd = ( nucSeqIndex )( std::rand( ) % ( pNucSeq->length( ) - iStart ) ) + iStart;

            bool bRev = ( std::rand( ) % 2 ) == 0;
            std::cout << iStart << " " << iEnd << ( bRev ? " true" : " false" ) << std::endl;

            if( bRev )
                pPack->vExtractSubsectionN( pPack->uiPositionToReverseStrand( iEnd ) + 1,
                                            pPack->uiPositionToReverseStrand( iStart ) + 1, xSeq2, false );
            else
                pPack->vExtractSubsectionN( iStart, iEnd, xSeq2, false );

            if( bRev )
            {
                xSeq2.vSwitchAllBasePairsToComplement( );
                xSeq2.vReverseAll( );
            } // if

            std::string sSeq3 = pNucSeq->fromTo( iStart, iEnd );

            std::cout << sSeq3 << std::endl;
            std::cout << xSeq2.toString( ) << std::endl;
            std::cout << sSeq3.compare( xSeq2.toString( ) ) << std::endl;
            std::cout << "----" << std::endl;

            assert( sSeq3.compare( xSeq2.toString( ) ) == 0 );
        } // for

    } // for

    return EXIT_SUCCESS;
} /// main function