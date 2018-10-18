#include "util/sqlite3.h"

/* Initialize the two transaction statements.
 */
void CppSQLiteDBExtended::vInititializeTransactionStatement( void )
{
    xStatementBeginTransaction = std::unique_ptr<CppSQLiteExtStatement>(
        new CppSQLiteExtStatement( *this, "BEGIN IMMEDIATE TRANSACTION" ) ); // precompiled transaction statement
    xStatementEndTransaction = std::unique_ptr<CppSQLiteExtStatement>(
        new CppSQLiteExtStatement( *this, "END TRANSACTION" ) ); // precompiled transaction statement
} // method

/* Creates the text for a SQL INSERT statement.
 * In the current version we automatically insert a NULL for the first column. (Primary Key)
 */
std::string sCreateSQLInsertStatementText(
    const char* pcTableName, // name of the table where we want to insert
    const unsigned int uiNumberOfArguemnts, // number of arguments for the table (-1, because we do not count the
                                            // primary key as first argument)
    bool bFirstColumnAsNULL // shall the first column automatically inserted as NULL
)
{
    std::string sInsertionStatement = "INSERT INTO ";
    sInsertionStatement.append( pcTableName ).append( bFirstColumnAsNULL ? " VALUES (NULL, " : " VALUES (" );

    /* Creation of comma separated list
     */
    for( unsigned int uiCounter = 0; uiCounter < uiNumberOfArguemnts; uiCounter++ )
    {
        sInsertionStatement.append( "@A" ).append( std::to_string( uiCounter ) );

        if( uiCounter < uiNumberOfArguemnts - 1 )
        {
            sInsertionStatement.append( ", " );
        } // if
    } // for

    sInsertionStatement.append( ")" );

    return sInsertionStatement;
} // function


/* Creates a fresh table within the database. */
void CppSQLiteDBExtended::vCreateTable(
    const std::vector<std::string>& xDatabaseColumns, // The names of the database columns
    const char* sTableName, // Table name
    bool bInsertColumnIdAsPrimaryKey // If we get a true here, then we insert a column id as primary key
)
{
    /* We drop the table in the case that it exists already. */
    execDML( std::string( "DROP TABLE IF EXISTS " ).append( sTableName ).c_str( ) );

    /* Generation of the table creation statement and its execution. */
    std::string sTableCreationStatement = "CREATE TABLE IF NOT EXISTS ";
    sTableCreationStatement.append( sTableName )
        .append( bInsertColumnIdAsPrimaryKey ? " (id INTEGER PRIMARY KEY, " : " (" );

    /* Create the sequence of column defintions in textual form */
    for( size_t i = 0; i < xDatabaseColumns.size( ) - 1; i++ )
        sTableCreationStatement.append( xDatabaseColumns[ i ] ).append( ", " );
    sTableCreationStatement.append( xDatabaseColumns.back( ) ).append( ")" );

    /* Execute the table statement */
    execDML( sTableCreationStatement.c_str( ) );
} // method


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

    CppSQLiteExtTable<int, std::string> table1(
        xTestDB, // the database where the table resides
        "table1", // name of the table in the database
        std::vector<std::string>{"number", "text"}, // column definitions of the table
        false // true == we build automatically a column for the primary key
    );
    {
        CppSQLiteExtImmediateTransactionContext xTransactionContext( xTestDB );
        std::cout << "START INSERT" << std::endl;
        table1.xInsertRow( 12, "abc " );
        table1.xInsertRow( 14, "gfd " );
        std::cout << "END INSERT" << std::endl;
    }
    table1.vDump( std::cout );
    table1.vForAllTableRowsUnpackedDo( []( int i, std::string s ) { std::cout << s << std::endl; } );

    CppSQLiteExtQueryStatement<int> statement ( xTestDB,
                                                 std::string( "SELECT number FROM " ).append( "table1" ).c_str( ) );
    statement.vExecuteAndForAllRowsUnpackedDo(  []( int i ) { std::cout << i << std::endl; });


} // function

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