#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#include "container/nucSeq.h"
#include <iostream>
#include <time.h>

#include <mysql_con.h>
#include <sql_api.h>
#include <db_con_pool.h>

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
            dDurtationComp += metaMeasureDuration( [&] { xCompressedNucSeq.compress( *pNucSeq ); } ).count( );
            NucSeq xNucSeqOut;
            dDurtationDeComp += metaMeasureDuration( [&] { xCompressedNucSeq.decompress( xNucSeqOut ); } ).count( );

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
        json{{SCHEMA, {{NAME, "nucseq_test"}, {FLAGS, {DROP_ON_CLOSURE}}}},
             {CONNECTION, {{HOSTNAME, "localhost"}, {USER, "root"}, {PASSWORD, "admin"}, {PORT, 0}}}} );

    json xCompNucSeqTableDef =
        json{{TABLE_NAME, "nucseq_table"}, {TABLE_COLUMNS, {{{COLUMN_NAME, "nucseqs"}}, {{COLUMN_NAME, "ins_id"}}}}};

    SQLTableWithAutoPriKey<SQLDB<MySQLConDB>, std::shared_ptr<CompressedNucSeq>, size_t> xTestTable(
        pDBCon,
        xCompNucSeqTableDef ); //
    xTestTable.deleteAllRows( ); // clear the table
    //using CompNucSeqSharedPtr = std::shared_ptr<CompressedNucSeq>;

    std::vector<std::shared_ptr<NucSeq>> vCompNucSeqs;
    for( size_t uiItr = 0; uiItr < 1000; uiItr++ )
    {
        auto pNucSeq = randomNucSeq( 250, 50, 10 );
        vCompNucSeqs.emplace_back( pNucSeq );
    } // for

    {
        auto pTrxnGuard = pDBCon->uniqueGuardedTrxn( );
        auto xBulkInserter = xTestTable.template getBulkInserter<300>( );
        metaMeasureAndLogDuration<true>( "Time required for table insertion:", //
                                         ( [&] {
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

void bulkInsertTest()
{
    double dLibIncrTotalTime = 0;
    double dAutoTotalTime = 0;
    for( int a = 0; a < 10; a++ )
        doNoExcept( [&] {
            json xTestTableLibIncrDef =
                json{{TABLE_NAME, "TEST_TABLE_LIB_INCR"},
                     {TABLE_COLUMNS, {{{COLUMN_NAME, "int_col_1"}}, {{COLUMN_NAME, "read"}}}},
                     /* { CPP_EXTRA, "DROP ON DESTRUCTION" } */};

            json xTestTableForeKeyDef =
                json{{TABLE_NAME, "TEST_TABLE_FOREIGN_KEY"},
                     {TABLE_COLUMNS, {{{COLUMN_NAME, "int_col_1"}}, {{COLUMN_NAME, "int_col_2"}}}},
                     /* { CPP_EXTRA, "DROP ON DESTRUCTION" } */};

            json xTestTableAutoDef = json{{TABLE_NAME, "TEST_TABLE_AUTO"},
                                          {TABLE_COLUMNS, {{{COLUMN_NAME, "int_col_1"}}, {{COLUMN_NAME, "read"}}}},
                                          /* { CPP_EXTRA, "DROP ON DESTRUCTION" } */};

            int numValues = 1000;

            std::vector<std::shared_ptr<NucSeq>> vCompNucSeqs1;
            std::vector<std::shared_ptr<NucSeq>> vCompNucSeqs2;
            for( size_t uiItr = 0; uiItr < (size_t)numValues; uiItr++ )
            {
                vCompNucSeqs1.emplace_back( randomNucSeq( 250, 0, 0 ) );
                vCompNucSeqs2.emplace_back( randomNucSeq( 250, 0, 0 ) );
            } // for

            double dLibIncrMaxTime = 0;
            double dAutoMaxTime = 0;
            std::vector<std::future<void>> vFutures;
            std::mutex xLock;
            {
                {
                    SQLDBConPool<MySQLConDB> xDBPool( 10, json{{SCHEMA, "test_pri_key"}} );
                    for( int i = 0; i < 10; i++ )
                        // type behind auto: std::shared_ptr<SQLDBConPool<MySQLConDB>::PooledSQLDBCon>
                        vFutures.push_back( xDBPool.enqueue( [&]( auto pDBCon ) {
                            doNoExcept(
                                [&] {
                                    SQLTableWithLibIncrPriKey<PooledSQLDBCon<MySQLConDB>,
                                                              int,
                                                              std::shared_ptr<CompressedNucSeq>>
                                        xPoolTable( pDBCon, xTestTableLibIncrDef );
                                    SQLTableWithLibIncrPriKey<PooledSQLDBCon<MySQLConDB>, PriKeyDefaultType,
                                                              PriKeyDefaultType>
                                        xPoolTableForeignKey( pDBCon, xTestTableForeKeyDef );
                                    auto pTrxnGuard = pDBCon->uniqueGuardedTrxn( );
                                    std::map<int, int> xVerificationMap;
                                    double dTime =
                                        metaMeasureAndLogDuration<true>( "BulkInserter required time:", [&]( ) {
                                            auto xBulkInserter = xPoolTable.template getBulkInserter<500>( );
                                            auto xBulkInserter2 = xPoolTableForeignKey.template getBulkInserter<500>( );
                                            std::lock_guard<std::mutex> xGuard( xLock );
                                            for( int i = 0; i < numValues; i++ )
                                            {
                                                auto uiPriKey1 = xBulkInserter->insert(
                                                    i * 2, makeSharedCompNucSeq( *vCompNucSeqs1[ i ] ) );
                                                auto uiPriKey2 = xBulkInserter->insert(
                                                    i * 2, makeSharedCompNucSeq( *vCompNucSeqs2[ i ] ) );
                                                //xBulkInserter2->insert( uiPriKey1, uiPriKey2 );
                                                // xVerificationMap[ uiPriKey ] = i * 2;
                                            }
                                            std::cout << "Finished inserting via BulkInserter LibIncr .... "
                                                      << std::endl;
                                        } );
                                    pDBCon->doPoolSafe( [&] { dLibIncrMaxTime = std::max( dLibIncrMaxTime, dTime ); } );

                                    // Verify the xVerificationMap
                                    // SQLQuery<PooledSQLDBCon<MySQLConDB>, int> xQuery(
                                    //     pDBCon, "SELECT int_col_1 FROM TEST_TABLE_LIB_INCR WHERE id = ?" );
                                    // for( auto rKeyValue : xVerificationMap )
                                    // {
                                    //     if( xQuery.scalar( rKeyValue.first ) != rKeyValue.second )
                                    //     {
                                    //         std::cout << "PANIC: Check failed" << std::endl;
                                    //         exit( 1 );
                                    //     } // if
                                    // } // for
                                },
                                "Problem during thread execution" );
                        } ) );
                }

                {
                    SQLDBConPool<MySQLConDB> xDBPool( 10, json{{SCHEMA, "test_pri_key"}} );
                    for( int i = 0; i < 0; i++ )
                        // type behind auto: std::shared_ptr<SQLDBConPool<MySQLConDB>::PooledSQLDBCon>
                        vFutures.push_back( xDBPool.enqueue( [&]( auto pDBCon ) {
                            doNoExcept(
                                [&] {
                                    SQLTableWithAutoPriKey<PooledSQLDBCon<MySQLConDB>,
                                                           int,
                                                           std::shared_ptr<CompressedNucSeq>>
                                        xPoolTable(
                                        pDBCon, xTestTableAutoDef );
                                    auto pTrxnGuard = pDBCon->uniqueGuardedTrxn( );

                                    // xPoolTable.insertNonSafe( 2, nullptr );
                                    // xPoolTable.insert( 2, 4 );
                                    double dTime =
                                        metaMeasureAndLogDuration<true>( "BulkInserter required time:", [&]( ) {
                                            auto xBulkInserter = xPoolTable.template getBulkInserter<500>( );
                                            std::lock_guard<std::mutex> xGuard( xLock );
                                            for( int i = 0; i < numValues; i++ )
                                                xBulkInserter->insert( i * 2,
                                                                       makeSharedCompNucSeq( *vCompNucSeqs1[ i ] ) );
                                            std::cout << "Finished inserting via BulkInserter Auto .... " << std::endl;
                                        } );
                                    pDBCon->doPoolSafe( [&] { dAutoMaxTime = std::max( dAutoMaxTime, dTime ); } );
                                },
                                "Problem during thread execution" );
                        } ) );
                }
            } // close the pool
            dLibIncrTotalTime += dLibIncrMaxTime;
            dAutoTotalTime += dAutoMaxTime;

            std::cout << "Finished one iteration .... " << std::endl;
        } ); // doNoExcept
    std::cout << "Libary increase: " << (int)( dLibIncrTotalTime / 10 ) << std::endl;
    std::cout << "MySQL increase: " << (int)( dAutoTotalTime / 10 ) << std::endl;
}

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

    // return doNoExcept( [&] { databaseTest( ); } );
    doNoExcept( [&] { bulkInsertTest( ); } );

    return EXIT_SUCCESS;
} /// main function
