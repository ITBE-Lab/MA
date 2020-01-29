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

#include <db_con_pool.h>
#include <common.h>
#include <threadPool.h>

using SQLStatement_ = SQLStatement<MySQLConDB>;
template <typename... Types> using SQLQuery_ = SQLQuery<MySQLConDB, Types...>;
template <typename... Types> using SQLTable_ = SQLTable<MySQLConDB, Types...>;


/* Meta function for measuring durations of function executions.
 * Using count() applied to the returned object we get the time in seconds.
 */
template <class FUNCTOR> std::chrono::duration<double> metaMeasureDuration( FUNCTOR&& f )
{ /* record start time
   */
    auto start = std::chrono::high_resolution_clock::now( );
    f( );

    /* record end time
     */
    auto end = std::chrono::high_resolution_clock::now( );
    return end - start;
} // meta function

template <class FUNCTOR>
void metaMeasureAndLogDuration( const std::string& sLogText, // additional logging text
                                FUNCTOR&& f // the functor called for measuring execution time
)
{
    auto xDuration = metaMeasureDuration( std::forward<FUNCTOR>( f ) );
    std::cout << sLogText << " required " << xDuration.count( ) * 1000 << " milliseconds." << std::endl;
} // meta

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
    print_tuple_impl( os, t, std::index_sequence_for<Args...>{ } );
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
    return a2t_impl( a, Indices{ } );
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
template <typename DBConnector> void checkDB( std::shared_ptr<SQLDB<DBConnector>> pMySQLDB, int jobnr )
{
    // std::cout << "Open Database" << std::endl;
    // std::shared_ptr<SQLDB<DBConnector>> pMySQLDB = std::make_shared<SQLDB<DBConnector>>( );

    {
        // json::array( { "WITHOUT ROWID" } )
        json xTestTableDef = {{TABLE_NAME, "TEST_TABLE" + std::to_string( jobnr )},
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
        // std::cout << std::setw( 2 ) << xTestTableDef << std::endl;
        // SQLTableWithAutoPriKey<DBConnector, int, double, SomeBlobType, std::string, uint32_t> xTestTable(
        //    pMySQLDB, xTestTableDef );
        SQLTableWithAutoPriKey<DBConnector, int, int, int, std::string, uint32_t> xTestTable( pMySQLDB, xTestTableDef );
        // xTestTable.deleteAllRows( );

        // std::array<std::tuple<std::nullptr_t, int, double, SomeBlobType, std::string, uint32_t>, 15> aArr;
        // for( int i = 0; i < 10; i++ )
        //     std::get<4>( aArr[ i ] ) = "This is a bit longer text that has to be written with every row";

		// LOAD DATA INFILE "C:\ProgramData\MySQL\MySQL Server 5.7\Uploads\0.csv" INTO TABLE test_table0;
        int numValues = 10000000 / 32;
        {
            auto pTrxnGuard = pMySQLDB->uniqueGuardedTrxn( );
            std::cout << "Make xBulkInserter" << std::endl;
			std::string text = "This"; // is a string of some size for insertion ...";
#if 1
            SomeBlobType blob;
            {
                metaMeasureAndLogDuration( "File Inserter required time:", [ & ]( ) {
                    auto xBulkInserter = xTestTable.template getFileBulkInserter<500>( std::to_string( jobnr ) );
                    for( int i = 0; i < numValues; i++ )
                        xBulkInserter->insert( nullptr, i, 4, i, text, std::numeric_limits<uint32_t>::max( ) );
                    std::cout << "Finished inserting .... " << std::endl;
                } );
            } // release of the bulk inserter
#endif

            //metaMeasureAndLogDuration( "Standard Bulk Inserter:", [ & ]( ) {
            //	auto xBulkNrmInserter = xTestTable.template getBulkInserter<500>();
            //    for( int i = 0; i < numValues; i++ )
            //        xBulkNrmInserter->insert( nullptr, i, 4, i, text, std::numeric_limits<uint32_t>::max( ) );
            //} );
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
        for( int32_t i = 1; i <= 10; i++ )
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
} // function


int main( int argc, char** argv )
{
    try
    {
        metaMeasureAndLogDuration( "Overall insertion time:", [ & ]( ) {
            
            SQLDBConPool<MySQLConDB> xDBPool( 32);
            for( int i = 0; i < 32; i++ )
                // type behind auto: std::shared_ptr<SQLDBConPool<MySQLConDB>::ConnectionManager>
                xDBPool.enqueue( []( auto pConMng ) {
                    std::cout << "In task: " << pConMng->getTaskId( ) << std::endl;
                    checkDB<MySQLConDB>( pConMng->pDBCon, (int)pConMng->getTaskId( ) );
                } );
        } );

#ifdef _MSC_VER
        int i;
        std::cin >> i;
#endif
        return 0;

        size_t uiPoolSize = 1;
        std::vector<std::shared_ptr<SQLDB<MySQLConDB>>> vec;
        for( unsigned int uiCount = 0; uiCount < uiPoolSize; uiCount++ )
        {
            vec.push_back( std::make_shared<SQLDB<MySQLConDB>>( ) );
        } // for

        {
            ThreadPool xPool( uiPoolSize ); // FIXME: Read this from parameters

            for( size_t uiJobId = 0; uiJobId < uiPoolSize; uiJobId++ )
                xPool.enqueue(
                    [ & ]( size_t, size_t uiJobId_ ) {
                        std::cout << "Start job Nr.: " << uiJobId_ << std::endl;
                        try
                        {
                            checkDB<MySQLConDB>( vec[ uiJobId_ ], (int)uiJobId_ );
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
        std::cout << "Database test failed. Reason:" << rxMySQLConExce.what( ) << std::endl;
        // return 1;
    } // catch
    catch( std::exception& rxException )
    {
        std::cout << "Database test failed. Reason: Runtime exception" << rxException.what( ) << std::endl;
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
