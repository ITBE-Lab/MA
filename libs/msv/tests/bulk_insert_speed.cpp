#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#include "ma/container/nucSeq.h"
#include <iostream>
#include <time.h>

#include <db_base.h>
#include <db_con_pool.h>
#include <sql_api.h>

#include "msv/container/sv_db/tables/pairedRead.h"
#include "msv/container/sv_db/tables/read.h"
#include "msv/container/sv_db/tables/sequencer.h"

using namespace libMSV;

std::shared_ptr<NucSeq> randomNucSeq( size_t uiLen, // length sequence
                                      size_t uiNMod = 200, // N start, statistically; 0 = no N at all
                                      size_t uiMaxSizeNseq = 150 ) // max size of N sequence

{
    auto pRet = std::make_shared<NucSeq>( );
    pRet->vReserveMemory( uiLen );
    int iNCounter = 0;
    for( size_t i = 0; i < uiLen; i++ )
    {
        pRet->sName = std::string( "randomly_generated_nucleotide_sequence_" ).append( std::to_string( i ) );
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

template <class T> std::string withCommas( T value )
{
    std::stringstream ss;
    ss.imbue( std::locale( "" ) );
    ss << std::fixed << value;
    return ss.str( );
} // method

void wait_for_everyone( std::mutex& xMutex,
                        size_t& uiThreadsReady,
                        std::condition_variable& xCondition,
                        size_t uiThreads )
{
    std::unique_lock<std::mutex> xLock( xMutex );
    uiThreadsReady++;
    if( uiThreadsReady < uiThreads + 1 )
        while( uiThreadsReady < uiThreads + 1 )
        {
            xCondition.wait( xLock );
        }
    else // the last thread that reaches the if will end up here
        xCondition.notify_all( );
} // method

int main( void )
{
    size_t uiThreads = std::thread::hardware_concurrency( );
    size_t uiNumValues = 19800;

    std::vector<std::shared_ptr<NucSeq>> vNucSeq1;
    std::vector<std::shared_ptr<NucSeq>> vNucSeq2;
    for( size_t uiItr = 0; uiItr < uiNumValues; uiItr++ )
    {
        vNucSeq1.emplace_back( randomNucSeq( 250, 0, 0 ) );
        vNucSeq2.emplace_back( randomNucSeq( 250, 0, 0 ) );
    } // for

    using duration = std::chrono::duration<double>;
    using time_point = std::chrono::time_point<std::chrono::steady_clock, duration>;


    std::condition_variable xCondition;
    std::mutex xMutex;
    size_t uiThreadsReady = 0;

    time_point xStart;
    {
        SQLDBConPool<DBConn> xDBPool(
            uiThreads, json{ { SCHEMA, { { NAME, "tmp_bulk_insert_speed" }, { FLAGS, { DROP_ON_CLOSURE } } } } } );

        std::vector<std::future<void>> vFutures;
        for( size_t i = 0; i < uiThreads; i++ )

            // type behind auto: std::shared_ptr<SQLDBConPool<DBConn>::PooledSQLDBCon>
            vFutures.push_back( xDBPool.enqueue( [ & ]( auto pDBCon ) {
                auto pSequencerTable = std::make_shared<SequencerTable<PooledSQLDBCon<DBConn>>>( pDBCon );
                auto pReadTable = std::make_shared<ReadTable<PooledSQLDBCon<DBConn>>>( pDBCon );
                auto pPairedReadTable = std::make_shared<PairedReadTable<PooledSQLDBCon<DBConn>>>( pDBCon, nullptr );
                auto pTrxnGuard = pDBCon->uniqueGuardedTrxn( );
                
                // The order of the BulkInserter is significant here.
                // We mist first destruct the xBulkInserter and then the xBulkInserter2
                auto xBulkInserter2 = pPairedReadTable->template getBulkInserter<1000>( );
                auto xBulkInserter = pReadTable->template getBulkInserter<500>();
                
                // wait till all threads are ready
                wait_for_everyone( xMutex, uiThreadsReady, xCondition, uiThreads );
                
                auto id = pSequencerTable->insert("Sequencer Run 1");
                for( size_t uiI = 0; uiI < uiNumValues; uiI++ )
                {
                    auto pKey1 =
                        xBulkInserter->insert( id, vNucSeq1[ uiI ]->sName, makeSharedCompNucSeq( *vNucSeq1[ uiI ] ) );
                    auto pKey2 =
                        xBulkInserter->insert( id, vNucSeq2[ uiI ]->sName, makeSharedCompNucSeq( *vNucSeq2[ uiI ] ) );
                    xBulkInserter2->insert( pKey1, pKey2 );
                } // for
            } ) );

        // wait till all threads reach this point
        wait_for_everyone( xMutex, uiThreadsReady, xCondition, uiThreads );
        xStart = std::chrono::steady_clock::now();

        // Check if something went wrong
        for( auto& rxFuture : vFutures )
        {
            try
            {
                rxFuture.get( );
            } // try
            catch (std::exception &rExcept)
            {
                std::cout << rExcept.what() << std::endl;
            } // catch
        } // for
    } // scope xDBPool -> wait for all threads to finish


    duration xTotalTime = std::chrono::steady_clock::now( ) - xStart;

    size_t uiTotalNumInserts = 3 * uiNumValues * uiThreads;
    size_t uiNumInsertedRows = (size_t)(uiTotalNumInserts / xTotalTime.count( ));
    std::cout << "inserted " << withCommas( uiNumInsertedRows ) << " rows per second (accumulated over " << uiThreads
              << " threads)." << std::endl;

    return EXIT_SUCCESS;
} /// main function
