#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#include "container/nucSeq.h"
#include <iostream>
#include <time.h>

#include <MySQL_con.h>
#include <common.h>

using namespace libMA;

std::shared_ptr<NucSeq> randomNucSeq( size_t uiLen, // length sequence
                                      size_t uiNMod = 200, // N start, statistically; 0 = no N at all
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

void compressedNucSeqTest( size_t uiLenStart, // start length sequence
                           size_t uiEndStart, // start length sequence
                           size_t uiNMod = 200, // statistically N start
                           size_t uiMaxSizeNseq = 150, // max size of N sequence
                           size_t uiTestRepeat = 5 ) // number of test repeat
{
    double dDurtationComp = 0;
    double dDurtationDeComp = 0;
    size_t uiCountMeasure = 0;
    std::cout << "Compression, decompression test ..." << std::endl;
    for( size_t iSize = uiLenStart; iSize < uiEndStart; iSize++ )
    {
        for( size_t iItr = 0; iItr < uiTestRepeat; iItr++ )
        {
            auto pNucSeq = randomNucSeq( iSize, uiNMod, uiMaxSizeNseq );
            CompressedNucSeq xCompressedNucSeq;
            uiCountMeasure++;
            dDurtationComp += metaMeasureDuration( [ & ] { xCompressedNucSeq.compress( *pNucSeq ); } ).count( );
            NucSeq xNucSeqOut;
            dDurtationDeComp += metaMeasureDuration( [ & ] { xCompressedNucSeq.decompress( xNucSeqOut ); } ).count( );

            if( pNucSeq->toString( ) != xNucSeqOut.toString( ) )
            {
                std::cout << "Discovered difference for sequence:" << xNucSeqOut.toString( ) << std::endl;
                exit( EXIT_FAILURE );
            } // if
            // DEBUG:
        } // for
    } // for
    std::cout << "Average compression time:" << ( dDurtationComp / uiCountMeasure ) * 1000 << std::endl;
    std::cout << "Average decompression time:" << ( dDurtationDeComp / uiCountMeasure ) * 1000 << std::endl;
} // function


void databaseTest( )
{
    std::shared_ptr<SQLDB<MySQLConDB>> pDBCon = std::make_shared<SQLDB<MySQLConDB>>(
        json{ { SCHEMA, "nucseq_test" },
              { CONNECTION, { { HOSTNAME, "localhost" }, { USER, "root" }, { PASSWORD, "admin" }, { PORT, 0 } } } } );

    json xCompNucSeqTableDef =
        json{ { TABLE_NAME, "nucseq_table" },
              { TABLE_COLUMNS, { { { COLUMN_NAME, "nucseqs" } }, { { COLUMN_NAME, "ins_id" } } } } };

    SQLTableWithAutoPriKey<SQLDB<MySQLConDB>, std::shared_ptr<CompressedNucSeq>, size_t> xTestTable(
        pDBCon,
        xCompNucSeqTableDef ); //
    xTestTable.deleteAllRows( ); // clear the table
    using CompNucSeqSharedPtr = std::shared_ptr<CompressedNucSeq>;

    std::vector<std::shared_ptr<NucSeq>> vCompNucSeqs;
    for( size_t uiItr = 0; uiItr < 10; uiItr++ )
    {
        auto pNucSeq = randomNucSeq( 25000000, 50, 10 );
        vCompNucSeqs.emplace_back( pNucSeq );
    } // for

    {
        auto pTrxnGuard = pDBCon->uniqueGuardedTrxn( );
        auto xBulkInserter = xTestTable.template getBulkInserter<300>( );
        metaMeasureAndLogDuration<true>( "Time required for table insertion:", //
                                         ( [ & ] {
                                             for( size_t uiItr = 0; uiItr < vCompNucSeqs.size( ); uiItr++ )
                                             {
                                                 xBulkInserter->insert(
                                                     nullptr, makeSharedCompNucSeq( *vCompNucSeqs[ uiItr ] ), uiItr );
                                             } // for
                                         } ) );
        std::cout << "Finished inserting .... " << std::endl;
    } // end transaction
#if 1
    for( size_t uiItr = 0; uiItr < vCompNucSeqs.size( ); uiItr++ )
    {
        SQLQuery<SQLDB<MySQLConDB>, std::shared_ptr<CompressedNucSeq>> xQuery(
            pDBCon, "SELECT nucseqs FROM nucseq_table WHERE ins_id = ?" );
        std::shared_ptr<NucSeq> pNucSeq = xQuery.execAndGetNthCell<0>( uiItr )->pUncomNucSeq;
        if( pNucSeq->toString( ) != vCompNucSeqs[ uiItr ]->toString( ) )
        {
            std::cout << "Discovered difference during retrieval:" << vCompNucSeqs[ uiItr ]->toString( ) << std::endl;
            exit( EXIT_FAILURE );
        } // if
    } // for
    std::cout << "Finished verification .... " << std::endl;
#endif
} // function

// Move Semantics
class A
{
  public:
    class B
    {
      public:
        B( const B& o )
        {
            std::cout << "Called Copy Constructor" << std::endl;
        }
        B( const B&& o )
        {
            std::cout << "Called Move Constructor" << std::endl;
        }

        B& operator=( B&& other )
        {
            std::cout << "Called Move Assignment" << std::endl;
            return *this;
        }

        B( )
        {}
        B( int i )
        {}
    };

  public:
    std::tuple<int, B> xTpl;
    B fun( )
    {
        return B( 2 );
    };

    void testTuple( int i, B&& rObj )
    {
        xTpl = std::tuple<int, B>( i, std::move( rObj ) );
    } // function
};

int main( void )
{
    // A::B xBObj;
    // A xAObj;
    // xAObj.testTuple( 2, std::move( xBObj ) );
    // std::cout << "Test end" << std::endl;
    // return 0;
    // 
    srand( static_cast<unsigned int>( time( NULL ) ) );
    // compressedNucSeqTest( 1, 100, 5, 5 ); // Test for short sequences with a mix of A,C,G,T,N
    // compressedNucSeqTest( 1, 100, 0, 0 ); // Test for short sequences with a mix of A,C,G,T only
    // compressedNucSeqTest( 1, 100, 1, 1000, 1 ); // Test for short N sequences
    // compressedNucSeqTest( 1000, 1500, 5, 50, 2 ); // Test for mid-sized sequences with a mix of A,C,G,T,N
    // compressedNucSeqTest( 1000, 1500, 0, 0, 2 ); // Test for mid-sized sequences with a mix of A,C,G,T only
    // compressedNucSeqTest( 1000, 1500, 1, 1000, 1 ); // Test for mid-sized N sequences
    // compressedNucSeqTest( 1000, 1500, 5, 100, 2 ); // Test for mid-sized sequences with a mix of A,C,G,T,N
    // compressedNucSeqTest( 1000, 1500, 0, 0, 2 ); // Test for mid-sized sequences with a mix of A,C,G,T only
    // compressedNucSeqTest( 1000, 1500, 1, 1000, 1 ); // Test for mid-sized N sequences
    // compressedNucSeqTest( 1000000, 1000065, 5, 500, 1 ); // Test for long sequences with a mix of A,C,G,T,N
    // compressedNucSeqTest( 1000000, 1000065, 0, 0, 1 ); // Test for long sequences with a mix of A,C,G,T only
    // compressedNucSeqTest( 1000000, 1000065, 1, 1000, 1 ); // Test for long N sequences
    // compressedNucSeqTest( 100000000, 100000001, 5, 500, 3 ); // Test for long sequences with a mix of A,C,G,T,N
    // compressedNucSeqTest( 100000000, 100000001, 0, 0, 3 ); // Test for long sequences with a mix of A,C,G,T only
    // compressedNucSeqTest( 100000000, 100000001, 1, 1000, 3 ); // Test for long N sequences

    return doNoExcept( [ & ] { databaseTest( ); } );

    return EXIT_SUCCESS;
} /// main function
