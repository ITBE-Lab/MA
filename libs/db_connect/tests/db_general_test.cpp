/*
 * src/test/examples/testlibpq.c
 *
 *
 * testlibpq.c
 *
 *      Test the C version of libpq, the PostgreSQL frontend library.
 */

// Then include the SQL common stuff that is currently below
// #define SQL_VERBOSE
#include <db_base.h>
#include <db_con_pool.h>
#include <sql_api.h>
#include <util/threadPool.h>

// #include "libpq-fe.h"
#include <stdio.h>
#include <stdlib.h>

// General definition of the DB connection
auto jDBConfig =
    json{ { SCHEMA, { { NAME, "db_test" } /* , { FLAGS, { DROP_ON_CLOSURE } } */ } },
// { TEMPORARY, true },
#ifdef USE_PG
          { CONNECTION, { { HOSTNAME, "localhost" }, { USER, "postgres" }, { PASSWORD, "admin" }, { PORT, 0 } } }
#endif
#ifdef USE_MSQL
          { CONNECTION, { { HOSTNAME, "localhost" }, { USER, "root" }, { PASSWORD, "admin" }, { PORT, 0 } } }
#endif
    };

template <typename DBConnector> void lowLevel( std::shared_ptr<DBConnector> pMySQLDB, size_t jobnr )
{
    // std::cout << "Open Database" << std::endl;
    // std::shared_ptr<SQLDB<DBConnector>> pMySQLDB = std::make_shared<SQLDB<DBConnector>>( );
    SQLStatement<DBConnector> xStmtDel( pMySQLDB, "DROP TABLE IF EXISTS TEST_1" );
    xStmtDel.exec( );

    SQLStatement<DBConnector> xStmt( pMySQLDB, "CREATE TABLE IF NOT EXISTS TEST_1 ( PersonID int8, Col2 int ) " );
    xStmt.exec( );

    SQLStatement<DBConnector> xInsertStmt( pMySQLDB, "INSERT INTO TEST_1(PersonID, Col2) VALUES($1, $2)" );
    for( int i = 0; i < 20; i++ )
        xInsertStmt.exec( (bool)i, i * 2 );
}

/* For checking with valgrind:
 * See: https://stackoverflow.com/questions/5134891/how-do-i-use-valgrind-to-find-memory-leaks
 * valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind-out.txt
 */
template <typename DBConnector> void checkDB( std::shared_ptr<DBConnector> pMySQLDB, size_t jobnr )
{


    {
        // json::array( { "WITHOUT ROWID" } )
        json xTestTableDef = json{ { TABLE_NAME, "TEST_TABLE" + std::to_string( jobnr ) },
                                   { TABLE_COLUMNS,
                                     { { { COLUMN_NAME, "Column_1" } /*, { CONSTRAINTS, "UNIQUE" } */ },
                                       { { COLUMN_NAME, "Column_2" } },
                                       { { COLUMN_NAME, "BlobColum" } },
                                       { { COLUMN_NAME, "TextColumn" } },
                                       { { COLUMN_NAME, "Col_uint32_t" } } } },
                                   { SQLITE_EXTRA, { "WITHOUT ROWID" } },
                                   { CPP_EXTRA, { "DROP ON DESTRUCTION" } }
                                   /* { SQL_EXTRA, { "INSERT NULL ON", 3 } } */ }; // ,
        // {CPP_EXTRA, "DROP ON DESTRUCTION"}};
        // std::cout << std::setw( 2 ) << xTestTableLibIncrDef << std::endl;
        SQLTableWithLibIncrPriKey<DBConnector, int, int, int, std::string, uint32_t> xTestTable( pMySQLDB,
                                                                                                 xTestTableDef );
        // Test next: SQLTableWithAutoPriKey
        // xTestTable.deleteAllRows( );

        int numValues = 100;
        {
            auto pTrxnGuard = pMySQLDB->uniqueGuardedTrxn( );
            std::cout << "Make xBulkInserter" << std::endl;
            std::string text = "This"; // is a string of some size for insertion ...";
            // text.resize( 20000 );

            // SomeBlobType blob;
            {
                metaMeasureAndLogDuration<true>( "BulkInserter required time:", [ & ]( ) {
                    auto xBulkInserter = xTestTable.template getBulkInserter<1000>( );
                    for( int i = 0; i < numValues; i++ )
                        xBulkInserter->insert( i, 4, i, text, std::numeric_limits<uint32_t>::max( ) );
                    std::cout << "Finished inserting .... " << std::endl;
                } );
            } // release of the bulk inserter

            //-- SQLQuery<DBConnector, int> xQueryScalar(pMySQLDB, "SELECT COUNT(*) FROM TEST_TABLE");
            //-- std::cout << "Count:" << xQueryScalar.scalar() << std::endl;
            // metaMeasureAndLogDuration( "Standard Bulk Inserter:", [ & ]( ) {
            //     auto xBulkNrmInserter = xTestTable.template getBulkInserter<300>( );
            //     for( int i = 0; i < numValues; i++ )
            //         xBulkNrmInserter->insert( nullptr, i, 4, i, text, std::numeric_limits<uint32_t>::max( ) );
            // } );
        } // scope transaction

        // json::array( { "WITHOUT ROWID" } )
        json xTestTableAutoIncrDef = json{ { TABLE_NAME, "TEST_TABLE_AUOINCR" + std::to_string( jobnr ) },
                                           { TABLE_COLUMNS,
                                             {
                                                 { { COLUMN_NAME, "Column1" } /*, { CONSTRAINTS, "UNIQUE" } */ },
                                                 { { COLUMN_NAME, "Column2" } /*, { CONSTRAINTS, "UNIQUE" } */ },
                                             } },
                                           { SQLITE_EXTRA, { "WITHOUT ROWID" } }
                                           // , { CPP_EXTRA, { "DROP ON DESTRUCTION" } }
                                           /* { SQL_EXTRA, { "INSERT NULL ON", 3 } } */ }; // ,
        // {CPP_EXTRA, "DROP ON DESTRUCTION"}};
        // std::cout << std::setw( 2 ) << xTestTableLibIncrDef << std::endl;
        SQLTableWithAutoPriKey<DBConnector, double, std::string> xTestTableAutoIncr( pMySQLDB, xTestTableAutoIncrDef );

        // SomeBlobType blob;
        {
            auto pTrxnGuard = pMySQLDB->uniqueGuardedTrxn( );
            metaMeasureAndLogDuration<true>( "BulkInserter required time:", [ & ]( ) {
                auto xBulkInserter = xTestTableAutoIncr.template getBulkInserter<10>( );
                for( int i = 0; i < numValues; i++ )
                    xBulkInserter->insert( 3.5, "Hello World" );
                std::cout << "Finished inserting .... " << std::endl;
            } );
        } // release of the bulk inserter

        SQLQuery<DBConnector, double, std::string> xQuery( pMySQLDB, "SELECT Column1, Column2 FROM TEST_TABLE_AUOINCR" +
                                                                         std::to_string( jobnr ) + " WHERE id=?" );
        xQuery.execAndForAll(
            [ & ]( double iCell, const std::string& rString ) { std::cout << iCell << " " << rString << std::endl; }, 2
            /* , std::string("Column1") */ );

        std::cout << "Do ImmediateQueryTmpl Test:" << std::endl;
        ImmediateQueryTmpl<std::shared_ptr<DBConnector>> xImmediateQuery(
            pMySQLDB, "SELECT Column1, Column2 FROM TEST_TABLE_AUOINCR" + std::to_string( jobnr ) );
        auto xTbl = xImmediateQuery.exec( );
        xTbl.print( );

        QueryExplainer<std::shared_ptr<DBConnector>> xQueryExplainer( pMySQLDB );
        xQueryExplainer.explain( "SELECT Column1, Column2 FROM TEST_TABLE_AUOINCR" + std::to_string( jobnr ) );
    };
#if 0
    }

        pMySQLDB
            ->addTable<int, int, int, std::string, uint32_t>(
                json{ { TABLE_NAME, "TEST_TABLE" + std::to_string(jobnr) },
                      { TABLE_COLUMNS,
                        { { { COLUMN_NAME, "Column_1" } /*, { CONSTRAINTS, "UNIQUE" } */ },
                          { { COLUMN_NAME, "Column_2" } },
                          { { COLUMN_NAME, "BlobColum" } },
                          { { COLUMN_NAME, "TextColumn" } },
                          { { COLUMN_NAME, "Col_uint32_t" } } } },
                      { SQLITE_EXTRA, { "WITHOUT ROWID" } } })
            .addTable<int, int, int, std::string, uint32_t>(
                json{ { TABLE_NAME, "TEST_TABLE" + std::to_string(jobnr) },
                      { TABLE_COLUMNS,
                        { { { COLUMN_NAME, "Column_1" } /*, { CONSTRAINTS, "UNIQUE" } */ },
                          { { COLUMN_NAME, "Column_2" } },
                          { { COLUMN_NAME, "BlobColum" } },
                          { { COLUMN_NAME, "TextColumn" } },
                          { { COLUMN_NAME, "Col_uint32_t" } } } },
                      { SQLITE_EXTRA, { "WITHOUT ROWID" } } });


        template <typename DBCon>
        using ReadTableType = SQLTableWithAutoPriKey<DBCon,
            int64_t, // sequencer id (foreign key)
            std::string, // read name
            NucSeqSql // read sequence
        >;
        const json jReadTableDef = {
            { TABLE_NAME, "read_table" },
            { TABLE_COLUMNS,
              { { { COLUMN_NAME, "sequencer_id" } }, { { COLUMN_NAME, "name" } }, { { COLUMN_NAME, "sequence" } } } },
            { FOREIGN_KEY, { { COLUMN_NAME, "sequencer_id" }, { REFERENCES, "sequencer_table(id)" } } } };


        auto xTblDescr(
            SQLTableDescr<DBCon,
            int64_t, // sequencer id (foreign key)
            std::string, // read name
            NucSeqSql> // read sequence
            ({ { TABLE_NAME, "read_table" },
                { TABLE_COLUMNS,
                  { { { COLUMN_NAME, "sequencer_id" } },
                    { { COLUMN_NAME, "name" } },
                    { { COLUMN_NAME, "sequence" } } } },
                { FOREIGN_KEY, { { COLUMN_NAME, "sequencer_id" }, { REFERENCES, "sequencer_table(id)" } } } }));

        xTblDescr::TableType


        // std::array<std::tuple<std::nullptr_t, int, double, SomeBlobType, std::string, uint32_t>, 15> aArr;
        // for( int i = 0; i < 10; i++ )
        //     std::get<4>( aArr[ i ] ) = "This is a bit longer text that has to be written with every row";

        // LOAD DATA INFILE "C:\ProgramData\MySQL\MySQL Server 5.7\Uploads\0.csv" INTO TABLE test_table0;
        int numValues = 1000;
        {
            auto pTrxnGuard = pMySQLDB->uniqueGuardedTrxn( );
            std::cout << "Make xBulkInserter" << std::endl;
            std::string text = "This"; // is a string of some size for insertion ...";
            // text.resize( 20000 );

            SomeBlobType blob;
            {
                // metaMeasureAndLogDuration<true>( "FileBulkInserter required time:", [ & ]( ) {
                //     auto xBulkInserter = xTestTable.template getFileBulkInserter<500>( );
                //     for( int i = 0; i < numValues; i++ )
                //         xBulkInserter->insert( nullptr, i, 4, i, text, std::numeric_limits<uint32_t>::max( ) );
                //     std::cout << "Finished inserting .... " << std::endl;
                // } );
            } // release of the bulk inserter

            //-- SQLQuery<DBConnector, int> xQueryScalar(pMySQLDB, "SELECT COUNT(*) FROM TEST_TABLE");
            //-- std::cout << "Count:" << xQueryScalar.scalar() << std::endl;
            // metaMeasureAndLogDuration( "Standard Bulk Inserter:", [ & ]( ) {
            //     auto xBulkNrmInserter = xTestTable.template getBulkInserter<300>( );
            //     for( int i = 0; i < numValues; i++ )
            //         xBulkNrmInserter->insert( nullptr, i, 4, i, text, std::numeric_limits<uint32_t>::max( ) );
            // } );
        } // scope transaction

 
        json xTestTable2Def = { { TABLE_NAME, "test_table_2" + std::to_string(jobnr) },
                                { TABLE_COLUMNS,
                                  { { { COLUMN_NAME, "int_column1" } /*, { CONSTRAINTS, "UNIQUE" } */ },
                                    { { COLUMN_NAME, "double_column" } },
                                    { { COLUMN_NAME, "int_column2" } } } } }; // ,
        SQLTable<DBConnector, uint64_t, double, uint32_t> xTable2(pMySQLDB, xTestTable2Def);
        {
            auto pTrxnGuard = pMySQLDB->uniqueGuardedTrxn();
            auto xBulkInserter = xTable2.template getBulkInserter<10>();

            metaMeasureAndLogDuration("Small Table: Time bulk insert:", [&]() {
                for (int i = 0; i < numValues; i++)
                    xBulkInserter->insert(1, 1, 1);
                });
        } // scope transaction

        json xTestTable3Def = { { TABLE_NAME, "test_table_3" + std::to_string(jobnr) },
                                { TABLE_COLUMNS,
                                  { { { COLUMN_NAME, "int_column1" } /*, { CONSTRAINTS, "UNIQUE" } */ },
                                    { { COLUMN_NAME, "double_column" } },
                                    { { COLUMN_NAME, "int_column2" } } } } }; // ,
        SQLTableWithAutoPriKey<DBConnector, uint64_t, double, uint32_t> xTable3(pMySQLDB, xTestTable3Def);
        {
            auto pTrxnGuard = pMySQLDB->uniqueGuardedTrxn();
            auto xBulkInserter = xTable3.template getBulkInserter<10>();

            metaMeasureAndLogDuration("Small Table: Time bulk insert:", [&]() {
                for (int i = 0; i < numValues; i++)
                    xBulkInserter->insert(nullptr, 1, 1, 1);
                });
        } // scope transaction


        return;
        // MySQLConDB::PreparedStmt obj( pMySQLDB, "INSERT" );
        // obj.bindAndExec(3.0, 4.0, "TExt");

        for( int32_t i = 1; i <= 2; i++ )
        {
            std::string text;
            SomeBlobType blob;
            for( int j = 0; j < i; j++ )
                text.append( "+" );
            auto val = xTestTable.insert( i + 10, i, 10, text, std::numeric_limits<uint32_t>::max( ) );
            std::cout << "Primary Key of prev insert: " << val << std::endl;
        } // for
        // std::cout << "Before NULL statement" << std::endl;
        // xTestTable.insertNonSafe( 1000, (double)5.0, nullptr, std::string( "Row with NULL" ) );
        // std::cout << "After NULL statement" << std::endl;

        // SQLStatement xTest( pMySQLDB, "SELECT * FROM TEST_TABLE" );
        xTestTable.dump( );

        // typename DBConnector::template PreparedQuery<int> xQuery( pMySQLDB, "SELECT Column_1 FROM TEST_TABLE" );

        SQLQuery_<int> xQuery( pMySQLDB, "SELECT Column_1 FROM TEST_TABLE" );
        xQuery.execAndForAll( [ & ]( int iCell ) { std::cout << iCell << std::endl; } );

        SQLQuery_<int> xQueryScalar( pMySQLDB, "SELECT COUNT(*) FROM TEST_TABLE" );
        std::cout << "Count:" << xQueryScalar.scalar( ) << std::endl;

        SQLQuery_<int64_t, int, double, SomeBlobType, std::string, uint32_t> xQuery2( pMySQLDB,
                                                                                      "SELECT * FROM TEST_TABLE" );
        std::cout << "Cell 0 : " << xQuery2.execAndGetNthCell<0>( ) << std::endl;
        std::cout << "Cell 1 : " << xQuery2.execAndGetNthCell<1>( ) << std::endl;
        std::cout << "Cell 2 : " << xQuery2.execAndGetNthCell<2>( ) << std::endl;
        std::cout << "Cell 3 : " << xQuery2.execAndGetNthCell<3>( ) << std::endl;
        std::cout << "Cell 4 : " << xQuery2.execAndGetNthCell<4>( ) << std::endl;
        std::cout << "Cell 5 OK : " << ( xQuery2.execAndGetNthCell<5>( ) == std::numeric_limits<uint32_t>::max( ) )
                  << std::endl;
        decltype( xTestTable )::CollTypeTranslation::dumpSQLColTypes( );

        pMySQLDB->execSQL( "CREATE INDEX text_index ON TEST_TABLE (Column_1)" );
        if( pMySQLDB->indexExists( "TEST_TABLE", "text_index" ) )
            std::cout << "INDEX text_index ON TEST_TABLE rediscovered!" << std::endl;

        xTestTable.addIndex( json{ /* { "INDEX_NAME", "sv_call_table_score_index_" } , */
                                   { "INDEX_COLUMNS", "Column_1" },
                                   { "WHERE", "Column_1 = 1" } } );
        xTestTable.addIndex( json{ /* { "INDEX_NAME", "sv_call_table_score_index_" } , */
                                   { "INDEX_COLUMNS", "Column_1" },
                                   { "WHERE", "Column_1 = 1" } } );

        std::cout << "Test table traversal via EOF" << std::endl;

        // xQuery.execAndBind( );
        // while( !xQuery.eof( ) )
        // {
        //     xQuery.next( );
        //     std::cout << std::get<0>( xQuery.get( ) ) << std::endl;
        // } // while
        //
        // std::cout << "Test table traversal via EOF 1" << std::endl;
        // xQuery.execAndBind( );
        // while( !xQuery.eof( ) )
        // {
        //     xQuery.next( );
        //     std::cout << std::get<0>( xQuery.get( ) ) << std::endl;
        // } // while
        //
        // std::cout << "Test table traversal via EOF 2" << std::endl;
        // xQuery.execAndBind( );
        // while( !xQuery.eof( ) )
        // {
        //     xQuery.next( );
        //     std::cout << std::get<0>( xQuery.get( ) ) << std::endl;
        // } // while

        // IMPORTANT TEST : executeAndStoreAllInVector, executeAndStoreInVector
        std::cout << xTestTable.makeInsertStmt( ) << std::endl;
        std::cout << xTestTable.makeInsertStmt( 3 ) << std::endl;
    }
#endif
} // function


void tableCreationTest1( std::shared_ptr<SQLDB<DBConn>> pDBConn, size_t jobnr )
{
    json xTestTableDef = json{ { TABLE_NAME, "TEST_TABLE" + std::to_string( jobnr ) },
                               { TABLE_COLUMNS,
                                 { { { COLUMN_NAME, "INT_COL_1" } /*, { CONSTRAINTS, "UNIQUE" } */ },
                                   { { COLUMN_NAME, "INT_COL_2" } },
                                   { { COLUMN_NAME, "INT_COL_3" } },
                                   { { COLUMN_NAME, "STRING_COL" } },
                                   { { COLUMN_NAME, "UINT32T_COL" } } } },
                               { CPP_EXTRA, { "DROP ON DESTRUCTION" } }
                               /* { SQL_EXTRA, { "INSERT NULL ON", 3 } } */ }; // ,
    SQLTableWithLibIncrPriKey<SQLDB<DBConn>, int, int, int, std::string, uint32_t> xTestTable( pDBConn, xTestTableDef );
} // test function

void tableCreationTest2(std::shared_ptr<SQLDB<DBConn>> pDBConn, size_t jobnr)
{
    json xTestTableDef = json{ { TABLE_NAME, "TEST_TABLE" },
                               { TABLE_COLUMNS,
                                 { { { COLUMN_NAME, "INT_COL_1" } },
                                   { { COLUMN_NAME, "INT_COL_2" } },
                                   { { COLUMN_NAME, "INT_COL_3" } },
                                   { { COLUMN_NAME, "STRING_COL" } },
                                   { { COLUMN_NAME, "UINT32T_COL" } } } } //,
                               //{ CPP_EXTRA, { "DROP ON DESTRUCTION" } }
    /* { SQL_EXTRA, { "INSERT NULL ON", 3 } } */ }; // ,
    SQLTableWithLibIncrPriKey<SQLDB<DBConn>, int, int, int, std::string, uint32_t> xTestTable(pDBConn, xTestTableDef);
} // test function

void tableCreationTest3(std::shared_ptr<SQLDB<DBConn>> pDBConn, size_t jobnr)
{
    json xTestTableDef = json{ { TABLE_NAME, "TEST_TABLE" },
                               { TABLE_COLUMNS,
                                 { { { COLUMN_NAME, "INT_COL_1" } },
                                   { { COLUMN_NAME, "INT_COL_2" } },
                                   { { COLUMN_NAME, "INT_COL_3" } },
                                   { { COLUMN_NAME, "STRING_COL" } },
                                   { { COLUMN_NAME, "UINT32T_COL" } } } },
                              { CPP_EXTRA, { "DROP ON DESTRUCTION" } }
    /* { SQL_EXTRA, { "INSERT NULL ON", 3 } } */ }; // ,
    SQLTableWithLibIncrPriKey<SQLDB<DBConn>, int, int, int, std::string, uint32_t> xTestTable(pDBConn, xTestTableDef);
} // test function

void textNullPtrInsertion(std::shared_ptr<SQLDB<DBConn>> pDBConn, size_t jobnr)
{
     json xTestTableDef = json{ { TABLE_NAME, "TEST_TABLE" },
                               { TABLE_COLUMNS,
                                 { { { COLUMN_NAME, "INT_COL_1" } },
                                   { { COLUMN_NAME, "INT_COL_2" } },
                                   { { COLUMN_NAME, "INT_COL_3" } } } } //,
                               //{ CPP_EXTRA, { "DROP ON DESTRUCTION" } }
    /* { SQL_EXTRA, { "INSERT NULL ON", 3 } } */ }; // ,
    SQLTableWithLibIncrPriKey<SQLDB<DBConn>, int, int, int> xTestTable(pDBConn, xTestTableDef);

    xTestTable.insertNonSafe(nullptr, 4, 6);
}

int test_in_pool( size_t uiPoolSize, const std::function<void( std::shared_ptr<SQLDB<DBConn>>, size_t )>& rfTestFun )
{
    try
    {
        std::vector<std::future<void>> vFutures;
        {
            SQLDBConPool<DBConn> xDBPool( uiPoolSize, jDBConfig );
            for( int i = 0; i < uiPoolSize; i++ )
                // type behind auto: std::shared_ptr<SQLDBConPool<MySQLConDB>::PooledSQLDBCon>
                vFutures.push_back( xDBPool.enqueue( [ &rfTestFun ]( auto pDBCon ) {
                    doNoExcept(
                        [ &rfTestFun, &pDBCon ] {
                            // pDBCon->doPoolSafe( [ & ] {
                            //     std::cout << "Thread:" << pDBCon->getTaskId( ) << " performs test" << std::endl;
                            // } );
                            // Database Test
                            rfTestFun( pDBCon, pDBCon->getTaskId( ) );
                        },
                        "Problem during thread execution" );
                } ) );
        } // close the pool

        // Get all futures for catching forwarded exceptions
        for( auto& rFurture : vFutures )
            doNoExcept( [ & ] { rFurture.get( ); } );
    } // try
    catch( std::runtime_error& e )
    {
        std::cout << "Test failed due to error:" << e.what( ) << std::endl;
        return 1;
    }
    return 0;
} // function

int main( int argc, char** argv )
{
    // Several connections create tables of different names
    test_in_pool( 1, textNullPtrInsertion );

    // Several connections create the same table
    //test_in_pool(20, tableCreationTest2);

    // Several connections create the same table
    //test_in_pool(1, tableCreationTest3);

    std::cout << "ALL WORK DONE ..." << std::endl;
    // #ifdef _MSC_VER
    //         int i;
    //         std::cin >> i;
    // #endif
    return 0;
#if 0
    size_t uiPoolSize = 1;
    std::vector<std::shared_ptr<SQLDB<DBConn>>> vec;
    for( unsigned int uiCount = 0; uiCount < uiPoolSize; uiCount++ )
    {
        vec.push_back( std::make_shared<SQLDB<DBConn>>(
            json{ { SCHEMA, { { NAME, "sv_db_2" }, { FLAGS, { DROP_ON_CLOSURE } } } } } ) );
    } // for

    {
        ThreadPool xPool( uiPoolSize ); // FIXME: Read this from parameters

        for( size_t uiJobId = 0; uiJobId < uiPoolSize; uiJobId++ )
            xPool.enqueue(
                [ & ]( size_t, size_t uiJobId_ ) {
                    std::cout << "Start job Nr.: " << uiJobId_ << std::endl; 
                    try
                    {
                        // checkDB( vec[ uiJobId_ ], uiJobId_ );
                        // Hmmm ...
                        // The concurrent construction of DBConnections seems to make trouble.
                        // std::shared_ptr<SQLDB<MySQLConDB>> pMySQLDB = std::make_shared<SQLDB<MySQLConDB>>();
                    }
                    catch( std::exception& e )
                    {
                        std::cout << "Job Nr.: " << uiJobId << " failed. Reason: " << e.what( ) << std::endl;
                    }
                },
                uiJobId );
    } // scope thread pool

} // try
#ifdef USE_PG
catch( PostgreSQLError& rxMySQLConExce )
{
    std::cout << "Database test failed. Reason:\n" << rxMySQLConExce.what( ) << std::endl;
    // return 1;
} // catch
#endif
catch( std::exception& rxException )
{
    std::cout << "Database test failed. Reason: Runtime exception:\n" << rxException.what( ) << std::endl;
    // return 1;
}
catch( ... )
{
    std::cout << "Database test failed. Reason: Some exception" << std::endl;
    //  return 1;
}

std::cout << "Database test finished successfully." << std::endl;

// std::array<std::tuple<int, CopyShower>, 40> aTestArray;
// tupleCatination( aTestArray );

// #ifdef _MSC_VER
//     int i;
//     std::cin >> i;
// #endif

return 0;
#endif
} // main
