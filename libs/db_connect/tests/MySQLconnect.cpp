/**
 * @brief MySQL connection test
 * @author Arne Kutzner, Markus Schmidt
 * @date Nov 2019
 */
#include <iostream>
// First include the connector
#include <MySQL_con.h>
// Then include the SQL common stuff that is currently below
// #define SQL_VERBOSE

#include <common.h>
#include <db_con_pool.h>
#include <threadPool.h>

using SQLStatement_ = SQLStatement<MySQLConDB>;
template <typename... Types> using SQLQuery_ = SQLQuery<MySQLConDB, Types...>;
template <typename... Types> using SQLTable_ = SQLTable<MySQLConDB, Types...>;

// Tuple printer for ostream.
// Found at: https://en.cppreference.com/w/cpp/utility/integer_sequence
template <class Ch, class Tr, class Tuple, std::size_t... Is>
void print_tuple_impl( std::basic_ostream<Ch, Tr>& os, const Tuple& t, std::index_sequence<Is...> )
{
    ( ( os << ( Is == 0 ? "" : ", " ) << std::get<Is>( t ) ), ... );
} // meta

template <class Ch, class Tr, class... Args>
auto& operator<<( std::basic_ostream<Ch, Tr>& os, const std::tuple<Args...>& t )
{
    os << "(";
    print_tuple_impl( os, t, std::index_sequence_for<Args...>{} );
    return os << ")";
} // meta

// Convert array into a tuple
// From: https://en.cppreference.com/w/cpp/utility/integer_sequence
template <typename Array, std::size_t... I> auto a2t_impl( const Array& a, std::index_sequence<I...> )
{
    return std::make_tuple( a[ I ]... );
}

template <typename T, std::size_t N, typename Indices = std::make_index_sequence<N>>
auto a2t( const std::array<T, N>& a )
{
    return a2t_impl( a, Indices{} );
}


template <typename... ArgTypes> void printArgPack( ArgTypes&&... args )
{
    std::stringstream xInfStream;
    ( ( xInfStream << typeid( args ).name( ) << ":" << args << "  " ), ... );
    std::cout << "Pack:" << xInfStream.str( ) << std::endl;
}

class CopyShower
{
  public:
    CopyShower( )
    {}
    CopyShower( const CopyShower& a )
    {
        std::cout << "MADE COPY" << std::endl;
    }
    CopyShower( const CopyShower&& a )
    {
        std::cout << "MADE MOVE" << std::endl;
    }
}; // class
template <class Ch, class Tr> auto& operator<<( std::basic_ostream<Ch, Tr>& os, const CopyShower& o )
{
    os << "+++" << std::endl;
    return os;
} // overloaded operator

template <size_t SIZE, class... TupleTypes>
void tupleCatination( const std::array<std::tuple<TupleTypes...>, SIZE>& aArr )
{
    auto tCatenated = tup_cat_arr( aArr );
    // using tplType = decltype( tup_cat_arr( aArr ) );
    std::cout << typeid( tCatenated ).name( ) << std::endl;
    std::cout << "START *********" << std::endl;
    STD_APPLY(
        [&]( auto&... args ) { // We forward the elements by reference for avoiding copies.
            printArgPack( args... ); // Here should occur the binding call
            // std::cout << "CONT" << std::endl;
        },
        tCatenated );
    std::cout << "END *********" << std::endl;
    std::cout << tCatenated << std::endl;
} // function

class SomeBlobType
{
  public:
    std::vector<char> vBuffer;

    SomeBlobType( ) : vBuffer( )
    {
        for( size_t uiCnt = 0; uiCnt < 20; uiCnt++ )
            vBuffer.push_back( char( uiCnt ) );
    } // constructor
}; // class

std::ostream& operator<<( std::ostream& os, const SomeBlobType& rBlob )
{
    os << "BLOB: ";
    for( char rVal : rBlob.vBuffer )
        os << (int)rVal << " ";
    return os;
} // operator

// Part1 : Specify the corresponding MySQL-type for your blob.
template <> inline std::string MySQLConDB::TypeTranslator::getSQLTypeName<SomeBlobType>( )
{
    return "LONGBLOB";
} // private function

// Part 2: Input arguments: Set the start of the blob (void *) and size of the blob.
template <> inline void MySQLConDB::StmtArg::set( const SomeBlobType& rxBlob )
{
    this->uiLength = static_cast<unsigned long>( rxBlob.vBuffer.size( ) );
    pMySQLBind->buffer_length = static_cast<unsigned long>( rxBlob.vBuffer.size( ) );
    pMySQLBind->buffer_type = MYSQL_TYPE_LONG_BLOB; // this type must be eqal to the type in Part 3.
    pMySQLBind->buffer = (void*)( &rxBlob.vBuffer[ 0 ] );
} // method

// Part 3: Code for supporting output of queries:
//         1. Via the third argument of init call set the MySQL datatype for your cell type.
//         2. Using storeVarSizeCel fetch the blob from the byte-buffer
template <> struct /* MySQLConDB:: */ RowCell<SomeBlobType> : public /* MySQLConDB:: */ RowCellBase<SomeBlobType>
{
    inline void init( MYSQL_BIND* pMySQLBind, SomeBlobType* pCellValue, size_t uiColNum )
    {
        RowCellBase<SomeBlobType>::init( pMySQLBind, pCellValue, MYSQL_TYPE_LONG_BLOB, uiColNum );
    } // method

    // Fetch the blob from the buffer.
    inline void storeVarSizeCell( )
    {
        // std::cout << "blob size: " << this->uiLength << std::endl;
        pCellValue->vBuffer.resize( this->uiLength );
        memcpy( (void*)&( pCellValue->vBuffer[ 0 ] ), this->pVarLenBuf.get( ), this->uiLength );
    } // method
}; // specialized class


/* For checking with valgrind:
 * See: https://stackoverflow.com/questions/5134891/how-do-i-use-valgrind-to-find-memory-leaks
 * valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind-out.txt
 */
template <typename DBConnector> void checkDB( std::shared_ptr<DBConnector> pMySQLDB, size_t jobnr )
{
    // std::cout << "Open Database" << std::endl;
    // std::shared_ptr<SQLDB<DBConnector>> pMySQLDB = std::make_shared<SQLDB<DBConnector>>( );

    {
        // json::array( { "WITHOUT ROWID" } )
        json xTestTableDef = json{{TABLE_NAME, "TEST_TABLE" + std::to_string( jobnr )},
                                  {TABLE_COLUMNS,
                                   {{{COLUMN_NAME, "Column_1"} /*, { CONSTRAINTS, "UNIQUE" } */},
                                    {{COLUMN_NAME, "Column_2"}},
                                    {{COLUMN_NAME, "BlobColum"}},
                                    {{COLUMN_NAME, "TextColumn"}},
                                    {{COLUMN_NAME, "Col_uint32_t"}}}},
                                  {SQLITE_EXTRA, {"WITHOUT ROWID"}} /*,
                              { CPP_EXTRA, { "DROP ON DESTRUCTION" } } */
                                  /* { SQL_EXTRA, { "INSERT NULL ON", 3 } } */}; // ,
        // {CPP_EXTRA, "DROP ON DESTRUCTION"}};
        // std::cout << std::setw( 2 ) << xTestTableLibIncrDef << std::endl;
        SQLTableWithAutoPriKey<DBConnector, int, int, int, std::string, uint32_t> xTestTable( pMySQLDB, xTestTableDef );
        // xTestTable.deleteAllRows( );

#if 0
        pMySQLDB
            ->addTable<int, int, int, std::string, uint32_t>(
                json{ { TABLE_NAME, "TEST_TABLE" + std::to_string( jobnr ) },
                      { TABLE_COLUMNS,
                        { { { COLUMN_NAME, "Column_1" } /*, { CONSTRAINTS, "UNIQUE" } */ },
                          { { COLUMN_NAME, "Column_2" } },
                          { { COLUMN_NAME, "BlobColum" } },
                          { { COLUMN_NAME, "TextColumn" } },
                          { { COLUMN_NAME, "Col_uint32_t" } } } },
                      { SQLITE_EXTRA, { "WITHOUT ROWID" } } } )
            .addTable<int, int, int, std::string, uint32_t>(
                json{ { TABLE_NAME, "TEST_TABLE" + std::to_string( jobnr ) },
                      { TABLE_COLUMNS,
                        { { { COLUMN_NAME, "Column_1" } /*, { CONSTRAINTS, "UNIQUE" } */ },
                          { { COLUMN_NAME, "Column_2" } },
                          { { COLUMN_NAME, "BlobColum" } },
                          { { COLUMN_NAME, "TextColumn" } },
                          { { COLUMN_NAME, "Col_uint32_t" } } } },
                      { SQLITE_EXTRA, { "WITHOUT ROWID" } } } );


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
            ( { { TABLE_NAME, "read_table" },
                { TABLE_COLUMNS,
                  { { { COLUMN_NAME, "sequencer_id" } },
                    { { COLUMN_NAME, "name" } },
                    { { COLUMN_NAME, "sequence" } } } },
                { FOREIGN_KEY, { { COLUMN_NAME, "sequencer_id" }, { REFERENCES, "sequencer_table(id)" } } } } ) );
		
		xTblDescr::TableType
#endif

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
#if 1
            SomeBlobType blob;
            {
                // metaMeasureAndLogDuration<true>( "FileBulkInserter required time:", [ & ]( ) {
                //     auto xBulkInserter = xTestTable.template getFileBulkInserter<500>( );
                //     for( int i = 0; i < numValues; i++ )
                //         xBulkInserter->insert( nullptr, i, 4, i, text, std::numeric_limits<uint32_t>::max( ) );
                //     std::cout << "Finished inserting .... " << std::endl;
                // } );
            } // release of the bulk inserter
#endif
            //-- SQLQuery<DBConnector, int> xQueryScalar(pMySQLDB, "SELECT COUNT(*) FROM TEST_TABLE");
            //-- std::cout << "Count:" << xQueryScalar.scalar() << std::endl;
            // metaMeasureAndLogDuration( "Standard Bulk Inserter:", [ & ]( ) {
            //     auto xBulkNrmInserter = xTestTable.template getBulkInserter<300>( );
            //     for( int i = 0; i < numValues; i++ )
            //         xBulkNrmInserter->insert( nullptr, i, 4, i, text, std::numeric_limits<uint32_t>::max( ) );
            // } );
        } // scope transaction

#if 0 
        json xTestTable2Def = { { TABLE_NAME, "test_table_2" + std::to_string( jobnr ) },
                                { TABLE_COLUMNS,
                                  { { { COLUMN_NAME, "int_column1" } /*, { CONSTRAINTS, "UNIQUE" } */ },
                                    { { COLUMN_NAME, "double_column" } },
                                    { { COLUMN_NAME, "int_column2" } } } } }; // ,
        SQLTable<DBConnector, uint64_t, double, uint32_t> xTable2( pMySQLDB, xTestTable2Def );
        {
            auto pTrxnGuard = pMySQLDB->uniqueGuardedTrxn( );
            auto xBulkInserter = xTable2.template getBulkInserter<10>( );

            metaMeasureAndLogDuration( "Small Table: Time bulk insert:", [ & ]( ) {
                for( int i = 0; i < numValues; i++ )
                    xBulkInserter->insert( 1, 1, 1 );
            } );
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
#endif

        return;
        // MySQLConDB::PreparedStmt obj( pMySQLDB, "INSERT" );
        // obj.bindAndExec(3.0, 4.0, "TExt");
#if 1
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
#endif
        // SQLStatement xTest( pMySQLDB, "SELECT * FROM TEST_TABLE" );
        xTestTable.dump( );

        // typename DBConnector::template PreparedQuery<int> xQuery( pMySQLDB, "SELECT Column_1 FROM TEST_TABLE" );

        SQLQuery_<int> xQuery( pMySQLDB, "SELECT Column_1 FROM TEST_TABLE" );
        xQuery.execAndForAll( [&]( int iCell ) { std::cout << iCell << std::endl; } );

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

        xTestTable.addIndex( json{/* { "INDEX_NAME", "sv_call_table_score_index_" } , */
                                  {"INDEX_COLUMNS", "Column_1"},
                                  {"WHERE", "Column_1 = 1"}} );
        xTestTable.addIndex( json{/* { "INDEX_NAME", "sv_call_table_score_index_" } , */
                                  {"INDEX_COLUMNS", "Column_1"},
                                  {"WHERE", "Column_1 = 1"}} );

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
} // function

/** @brief Executes the func so that exceptions are swallowed.
 *  @detail For destructors and threads so that they do not throw.
 */
template <typename F> inline void doSwallowingExcpt( F&& func )
{
    try
    {
        std::cout << "Do swallowing start..." << std::endl;
        func( );
        std::cout << "Do swallowing end..." << std::endl;
    } // try
    catch( std::exception& rxExcpt )
    {
        // Swallow the exception
        std::cout << std::string( "SQLDBConPool: Task/thread execution terminated abnormally. Details:\n" ) +
                         rxExcpt.what( )
                  << std::endl;
    } // catch
    catch( ... )
    {
        std::cout << "SQLDBConPool: Task/thread execution terminated abnormally for unknown reasons." << std::endl;
    } // catch
    std::cout << "Do swallowing terminated..." << std::endl;
} // method

void excptTest( )
{
#if 0
	std::cout << "Before do..." << std::endl;
	doSwallowingExcpt( [] {
        std::cout << "Before throw..." << std::endl;
        throw std::runtime_error( "A problem ..." );
		std::cout << "After throw..." << std::endl;
    } );
	std::cout << "Afer do..." << std::endl;
#endif
    SQLDBConPool<MySQLConDB> xDBPool( 1 );

    // type behind auto: std::shared_ptr<SQLDBConPool<MySQLConDB>::PooledSQLDBCon>
    auto xFuture = xDBPool.enqueue( []( auto pConMng ) { throw std::runtime_error( "Throw ..." ); } );

    doSwallowingExcpt( [&]( ) { xFuture.get( ); } ); // catch via future
} // method

/** @brief Test the primary key related part of the code */
int priKeyTest( )
{
    double dLibIncrTotalTime = 0;
    double dAutoTotalTime = 0;
    for( int a = 0; a < 10; a++ )
        doNoExcept( [&] {
            json xTestTableLibIncrDef =
                json{{TABLE_NAME, "TEST_TABLE_LIB_INCR"},
                     {TABLE_COLUMNS, {{{COLUMN_NAME, "int_col_1"}}, {{COLUMN_NAME, "int_col_2"}}}},
                     /* { CPP_EXTRA, "DROP ON DESTRUCTION" } */};

            json xTestTableAutoDef = json{{TABLE_NAME, "TEST_TABLE_AUTO"},
                                          {TABLE_COLUMNS, {{{COLUMN_NAME, "int_col_1"}}, {{COLUMN_NAME, "int_col_2"}}}},
                                          /* { CPP_EXTRA, "DROP ON DESTRUCTION" } */};

            int numValues = 10000;

            double dLibIncrMaxTime = 0;
            double dAutoMaxTime = 0;
            std::vector<std::future<void>> vFutures;
            // std::mutex xLock;
            {
                {
                    SQLDBConPool<MySQLConDB> xDBPool( 10, json{{SCHEMA, "test_pri_key"}} );
                    for( int i = 0; i < 20; i++ )
                        // type behind auto: std::shared_ptr<SQLDBConPool<MySQLConDB>::PooledSQLDBCon>
                        vFutures.push_back( xDBPool.enqueue( [&]( auto pDBCon ) {
                            doNoExcept(
                                [&] {
                                    SQLTableWithLibIncrPriKey<PooledSQLDBCon<MySQLConDB>, int, int> xPoolTable(
                                        pDBCon, xTestTableLibIncrDef );
                                    auto pTrxnGuard = pDBCon->uniqueGuardedTrxn( );
                                    std::map<int, int> xVerificationMap;
                                    double dTime =
                                        metaMeasureAndLogDuration<true>( "BulkInserter required time:", [&]( ) {
                                            auto xBulkInserter = xPoolTable.template getBulkInserter<500>(
                                                RESERVE_BATCH_OF_PRIKEY ); // RESERVE_BATCH_OF_PRIKEY );
                                            // std::lock_guard<std::mutex> xGuard( xLock );
                                            for( int i = 0; i < numValues; i++ )
                                            {
                                                auto uiPriKey = xBulkInserter->insert( i * 2, 10 );
                                                xVerificationMap[ uiPriKey ] = (int)( i * 2 );
                                            }
                                            std::cout << "Finished inserting via BulkInserter LibIncr .... "
                                                      << std::endl;
                                        } );
                                    pDBCon->doPoolSafe( [&] { dLibIncrMaxTime = std::max( dLibIncrMaxTime, dTime ); } );

                                    // Verify the xVerificationMap
                                    SQLQuery<PooledSQLDBCon<MySQLConDB>, int> xQuery(
                                        pDBCon, "SELECT int_col_1 FROM TEST_TABLE_LIB_INCR WHERE id = ?" );
                                    for( auto rKeyValue : xVerificationMap )
                                    {
                                        if( xQuery.scalar( rKeyValue.first ) != rKeyValue.second )
                                        {
                                            std::cout << "PANIC: Check failed" << std::endl;
                                            exit( 1 );
                                        } // if
                                    } // for
                                },
                                "Problem during thread execution" );
                        } ) );
                }
#if 0
                {
                    SQLDBConPool<MySQLConDB> xDBPool( 10, json{ { SCHEMA, "test_pri_key" } } );
                    for( int i = 0; i < 10; i++ )
                        // type behind auto: std::shared_ptr<SQLDBConPool<MySQLConDB>::PooledSQLDBCon>
                        vFutures.push_back( xDBPool.enqueue( [ & ]( auto pDBCon ) {
                            doNoExcept(
                                [ & ] {
                                    SQLTableWithAutoPriKey<PooledSQLDBCon<MySQLConDB>, int, int> xPoolTable(
                                        pDBCon, xTestTableAutoDef );
                                    auto pTrxnGuard = pDBCon->uniqueGuardedTrxn( );

                                    // xPoolTable.insertNonSafe( 2, nullptr );
                                    // xPoolTable.insert( 2, 4 );
                                    double dTime =
                                        metaMeasureAndLogDuration<true>( "BulkInserter required time:", [ & ]( ) {
                                            auto xBulkInserter = xPoolTable.template getBulkInserter<500>( );
                                            // std::lock_guard<std::mutex> xGuard( xLock );
                                            for( int i = 0; i < numValues; i++ )
                                                xBulkInserter->insert( i * 2, 10 );
                                            std::cout << "Finished inserting via BulkInserter Auto .... " << std::endl;
                                        } );
                                    pDBCon->doPoolSafe( [ & ] { dAutoMaxTime = std::max( dAutoMaxTime, dTime ); } );
                                },
                                "Problem during thread execution" );
                        } ) );
                }
#endif
            } // close the pool
            dLibIncrTotalTime += dLibIncrMaxTime;
            dAutoTotalTime += dAutoMaxTime;

            std::cout << "Finished one iteration .... " << std::endl;
        } ); // doNoExcept
    std::cout << "Libary increase: " << (int)( dLibIncrTotalTime / 10 ) << std::endl;
    std::cout << "MySQL increase: " << (int)( dAutoTotalTime / 10 ) << std::endl;
    // std::shared_ptr<SQLDB<MySQLConDB>> pDBCon =
    //     std::make_shared<SQLDB<MySQLConDB>>( json{ { SCHEMA, "test_pri_key" } } );
    // pDBCon->execSQL( "DROP TABLE TEST_TABLE" );

    return 0;
} // function (priKeyTest)


/** @brief Test for correct drop schema in the context of a connection pool */
int dropSchemaTest( )
{
    SQLDBConPool<MySQLConDB> xDBPool( 10, json{{SCHEMA, {{NAME, "test_schema_drop"}, {FLAGS, {DROP_ON_CLOSURE}}}}} );

    return 0;

} // function (dropSchemaTest)


/** @brief Test the listing of all schema */
int schemaListTest( )
{
    std::shared_ptr<SQLDB<MySQLConDB>> xDBCon( std::make_shared<SQLDB<MySQLConDB>>(
        json{{SCHEMA, {{NAME, "schema_list_test"}, {FLAGS, {DROP_ON_CLOSURE}}}}} ) );
    auto xDBInformer = SQLDBInformer<SQLDB<MySQLConDB>>( xDBCon );
    std::cout << "List of schema for the current connection" << std::endl;
    for( auto& sSchema : xDBInformer.getAllSchemas( ) )
        std::cout << sSchema << std::endl;

    return 0;
} // function (schemaListTest)


int mySQLStorageResultTest( )
{
    doNoExcept( [&] {
        std::cout << "Do mySQLStorageResultTest ..." << std::endl;
        std::shared_ptr<SQLDB<MySQLConDB>> pDBCon(
            std::make_shared<SQLDB<MySQLConDB>>( json{{SCHEMA, {{NAME, "storage_result_test"}}}} ) );

        json xTestTableWithStringsDef =
            json{{TABLE_NAME, "TEST_TABLE_WITH_STRINGS"}, {TABLE_COLUMNS, {{{COLUMN_NAME, "string_col"}}}}};

        SQLTableWithLibIncrPriKey<SQLDB<MySQLConDB>, std::string> xTestTableWithStrings( pDBCon,
                                                                                         xTestTableWithStringsDef );

        auto xBulkInserter = xTestTableWithStrings.template getBulkInserter<500>( RESERVE_BATCH_OF_PRIKEY );

        {
            auto pTrxnGuard = pDBCon->uniqueGuardedTrxn( );
            for( int i = 0; i < 10000000; i++ )
                xBulkInserter->insert( "This is a text that requires some buffer memory ..." );
        } // end of transaction

        // Verify the xVerificationMap
        std::cout << "Do execute query ..." << std::endl;
        SQLQuery<SQLDB<MySQLConDB>, std::string> xQuery(
            pDBCon,
            "SELECT string_col FROM TEST_TABLE_WITH_STRINGS",
            json{{MYSQL_EXTRA, {{FLAGS, {NO_RESULT_STORAGE_ON_CLIENT}}}}} );
        size_t uiCount = 0;
        xQuery.execAndForAll( [&]( auto iCell ) { uiCount++; } );

        std::cout << "Counter: " << uiCount << std::endl;
    } );

    return 0;
} // function (mySQLStorageResultTest)

int main( int argc, char** argv )
{
    priKeyTest( );
    return 0;

    try
    {
        // definition of database connection in json:
        auto jDBConfig = json{{SCHEMA, {{NAME, "sv_db"}, {FLAGS, {DROP_ON_CLOSURE}}}},
                              // { TEMPORARY, true },
                              {CONNECTION, {{HOSTNAME, "localhost"}, {USER, "root"}, {PASSWORD, "admin"}, {PORT, 0}}}};

        std::vector<std::future<void>> vFutures;
        {
            SQLDBConPool<MySQLConDB> xDBPool( 10, jDBConfig );
            for( int i = 0; i < 10 /* 20 */; i++ )
                // type behind auto: std::shared_ptr<SQLDBConPool<MySQLConDB>::PooledSQLDBCon>
                vFutures.push_back( xDBPool.enqueue( []( auto pDBCon ) {
                    doNoExcept(
                        [&] {
                            typedef decltype( *pDBCon ) Type;
                            // typedef typename SubType :: element_type Type;
                            // using Type = typename decltype( pDBCon )::element_type;
                            std::cout << "Job executed in task: " << pDBCon->getTaskId( ) << std::endl;
                            checkDB( pDBCon, pDBCon->getTaskId( ) );
                            pDBCon->doPoolSafe( [] { std::cout << "This print is pool safe ..." << std::endl; } );
                        },
                        "Problem during thread execution" );
                } ) );
        } // close the pool

        // Get all future exception safe
        for( auto& rFurture : vFutures )
            doNoExcept( [&] { rFurture.get( ); } );

        std::cout << "ALL WORK DONE ..." << std::endl;
#ifdef _MSC_VER
        int i;
        std::cin >> i;
#endif
        return 0;

        size_t uiPoolSize = 1;
        std::vector<std::shared_ptr<SQLDB<MySQLConDB>>> vec;
        for( unsigned int uiCount = 0; uiCount < uiPoolSize; uiCount++ )
        {
            vec.push_back( std::make_shared<SQLDB<MySQLConDB>>(
                json{{SCHEMA, {{NAME, "sv_db_2"}, {FLAGS, {DROP_ON_CLOSURE}}}}} ) );
        } // for

        {
            ThreadPool xPool( uiPoolSize ); // FIXME: Read this from parameters

            for( size_t uiJobId = 0; uiJobId < uiPoolSize; uiJobId++ )
                xPool.enqueue(
                    [&]( size_t, size_t uiJobId_ ) {
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
    catch( MySQLConException& rxMySQLConExce )
    {
        std::cout << "Database test failed. Reason:\n" << rxMySQLConExce.what( ) << std::endl;
        // return 1;
    } // catch
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

#ifdef _MSC_VER
    int i;
    std::cin >> i;
#endif

    return 0;
} // main
