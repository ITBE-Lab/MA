#include <db_sql.h>
// #include "container/nucSeq.h" // Required for test-function

#if 0
using namespace libMA;
/** @brief SQLite interface usage example function
 */
void doSqlite3Test( void )
{
    CppSQLiteDBExtended xTestDB( "./", // directory location of the database
                                 "test_db", // name of the database
                                 eCREATE_DB ); // DB create/open

    // /* We create the insert statement for the table alignment results.
    //  */
    // CppSQLiteExtInsertStatement<int, // 1. analysis_id INTEGER
    //                             const char*, // 2. sequence_name TEXT
    //                             >
    //     xSQLInsertAlignmentResults( rxDatabase, "alignment_results" );
    //
    // xTestDB.vCreateTable( std::vector<std::string>{"number INTEGER", // Foreign Key The gene_id of the Gene
    //                                                "text TEXT"}, // vector initialization
    //                       "table1" // sTableName
    // );

    NucSeqSql xSeq( std::make_shared<NucSeq>( "ACCGTG" ) );

    CppSQLiteExtTableWithAutomaticPrimaryKey<int64_t, NucSeqSql> table1(
        xTestDB, // the database where the table resides
        "table1", // name of the table in the database
        std::vector<std::string>{"number", "nucSeq"} // column definitions of the table
    );
    {
        CppSQLiteExtImmediateTransactionContext xTransactionContext( xTestDB );
        std::cout << "START INSERT" << std::endl;
        table1.xInsertRow( 12, xSeq );
        table1.xInsertRow( 14, xSeq );
        std::cout << "END INSERT" << std::endl;
    }
    table1.vDump( std::cout );
    table1.vForAllTableRowsUnpackedDo(
        []( int, int64_t i, NucSeqSql s ) { std::cout << s.pNucSeq->toString( ) << std::endl; } );

    CppSQLiteExtQueryStatement<int> statement( xTestDB, "SELECT number FROM table1" );
    statement.vExecuteAndForAllRowsUnpackedDo( []( int64_t i ) { std::cout << i << std::endl; } );


} // function
#endif

#ifdef SQLITE_MAIN
int main( void )
{
    try
    {
        doSqlite3Test( );
    }
    catch( const std::exception& e )
    {
        std::cout << e.what( ) << std::endl;
    }
    catch( ... )
    {
        std::cout << "Catched general exception" << std::endl;
    }
} // function
#endif