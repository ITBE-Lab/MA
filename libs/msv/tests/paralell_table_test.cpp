#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#include "ma/container/nucSeq.h"
#include <iostream>
#include <time.h>

#include <db_con_pool.h>
#include <mysql_con.h>
#include <sql_api.h>

#include "msv/container/sv_db/tables/pairedRead.h"
#include "msv/container/sv_db/tables/read.h"
#include "msv/container/sv_db/tables/sequencer.h"

using namespace libMSV;


template <typename DBCon>
using TestTableType = SQLTableWithLibIncrPriKey<DBCon,
                                                uint32_t, // pos
                                                WKBUint64Rectangle // rectangle (geometry)
                                                >;


template <typename DBCon> class TestTable : public TestTableType<DBCon>
{
  public:
    json jSvCallTableDef( size_t uiSuffix )
    {
        return json{{TABLE_NAME, "test_table_" + std::to_string( uiSuffix )},
                    {TABLE_COLUMNS, {{{COLUMN_NAME, "pos"}}, {{COLUMN_NAME, "rectangle"}, {CONSTRAINTS, "NOT NULL"}}}}};
    } // method

    TestTable( std::shared_ptr<DBCon> pDatabase, size_t uiSuffix )
        : TestTableType<DBCon>( pDatabase, // the database where the table resides
                                jSvCallTableDef( uiSuffix ) )

    {} // default constructor
};


int main( void )
{
    size_t uiNumValues = 6000000;

    size_t uiNumTables = 16;
    size_t uiMaxSizeRect = 10000;

    std::vector<std::vector<WKBUint64Rectangle>> vvRectangles;
    std::vector<std::vector<uint32_t>> vvPos;
    for( size_t uiI = 0; uiI < uiNumTables; uiI++ )
    {
        vvRectangles.emplace_back( );
        vvPos.emplace_back( );
        for( size_t uiJ = 0; uiJ < uiNumValues / uiNumTables; uiJ++ )
        {
            auto xRectangle = WKBUint64Rectangle( geom::Rectangle<nucSeqIndex>( std::rand( ) % uiMaxSizeRect,
                                                                                std::rand( ) % uiMaxSizeRect,
                                                                                std::rand( ) % uiMaxSizeRect,
                                                                                std::rand( ) % uiMaxSizeRect ) );
            vvRectangles.back( ).push_back( xRectangle );
            vvPos.back( ).push_back( std::rand( ) % uiMaxSizeRect );
        } // for
    } // for

    using duration = std::chrono::duration<double>;
    using time_point = std::chrono::time_point<std::chrono::steady_clock, duration>;


    time_point xStart;

    SQLDBConPool<MySQLConDB> xDBPool( uiNumTables,
                                      json{{SCHEMA, {{NAME, "tmp_parallel_table_test"}, {FLAGS, {DROP_ON_CLOSURE}}}}} );
    std::vector<std::future<void>> vFutures;

    xStart = std::chrono::steady_clock::now( );
    for( size_t i = 0; i < uiNumTables; i++ )
        // type behind auto: std::shared_ptr<SQLDBConPool<MySQLConDB>::PooledSQLDBCon>
        vFutures.push_back( xDBPool.enqueue(
            [&]( auto pDBCon, size_t iThread ) {
                auto pTestTable = std::make_shared<TestTable<PooledSQLDBCon<MySQLConDB>>>( pDBCon, iThread );

                auto pTrxnGuard = pDBCon->uniqueGuardedTrxn( );

                auto xBulkInserter = pTestTable->template getBulkInserter<500>( );

                for( size_t uiI = 0; uiI < uiNumValues / uiNumTables; uiI++ )
                {
                    xBulkInserter->insert( vvPos[ iThread ][ uiI ], vvRectangles[ iThread ][ uiI ] );
                } // for
            },
            i ) );


    // Get all future exception safe
    for( auto& rFurture : vFutures )
        doNoExcept( [&] { rFurture.get( ); } );

    duration xTotalTime = std::chrono::steady_clock::now( ) - xStart;

    std::cout << "inserted " << uiNumValues << " rows in " << xTotalTime.count( ) << " seconds using " << uiNumTables
              << " tables." << std::endl;


    vFutures.clear( );

    xStart = std::chrono::steady_clock::now( );
    for( size_t i = 0; i < uiNumTables; i++ )
        // type behind auto: std::shared_ptr<SQLDBConPool<MySQLConDB>::PooledSQLDBCon>
        vFutures.push_back( xDBPool.enqueue(
            [&]( auto pDBCon, size_t iThread ) {
                auto pTestTable = std::make_shared<TestTable<PooledSQLDBCon<MySQLConDB>>>( pDBCon, iThread );

                std::cout << "triggering index construction " << iThread << std::endl;
                pTestTable->addIndex(
                    json{{INDEX_NAME, "rectangle"}, {INDEX_COLUMNS, "rectangle"}, {INDEX_TYPE, "SPATIAL"}}, false );
                //pTestTable->addIndex(
                //    json{{INDEX_NAME, "pos_index" + std::to_string( iThread )}, {INDEX_COLUMNS, "pos"}} );
                std::cout << "finished index construction " << iThread << std::endl;
            },
            i ) );


    // Get all future exception safe
    for( auto& rFurture : vFutures )
        doNoExcept( [&] { rFurture.get( ); } );

    xTotalTime = std::chrono::steady_clock::now( ) - xStart;
    std::cout << "built " << uiNumTables << " indices in " << xTotalTime.count( ) << " seconds" << std::endl;

    //int c;
    //std::cin >> c;
    return EXIT_SUCCESS;
} /// main function
