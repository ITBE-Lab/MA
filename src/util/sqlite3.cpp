#include "util/sqlite3.h"
#include "container/nucSeq.h"

// list of valied types:
template <> std::string getSQLTypeName<std::string>( )
{
    return "TEXT";
}
// numeric:
// full number:
// sqLITE maps all numbers to INTEGER anyways?
// signed
template <> std::string getSQLTypeName<int8_t>( )
{
    return "INTEGER";
}
template <> std::string getSQLTypeName<int16_t>( )
{
    return "INTEGER";
}
template <> std::string getSQLTypeName<int32_t>( )
{
    return "INTEGER";
}
template <> std::string getSQLTypeName<int64_t>( )
{
    return "INTEGER";
}

// unsigned
template <> std::string getSQLTypeName<uint8_t>( )
{
    return "INTEGER";
}
template <> std::string getSQLTypeName<uint16_t>( )
{
    return "INTEGER";
}
template <> std::string getSQLTypeName<uint32_t>( )
{
    return "INTEGER";
}

// template <> std::string SQLTypeName<uint64_t>::name = "INTEGER"; <- this is too large for sqlite3

// floating point:
template <> std::string getSQLTypeName<double>( )
{
    return "DOUBLE";
}
template <> std::string getSQLTypeName<float>( )
{
    return "FLOAT";
}

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
    std::cout << sTableCreationStatement << std::endl;
    execDML( sTableCreationStatement.c_str( ) );
} // method

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

    NucSeqSql xSeq;
    xSeq.pNucSeq = std::make_shared<NucSeq>( "ACCGTG" );

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
    table1.vForAllTableRowsUnpackedDo( []( int, int64_t i, NucSeqSql s ) { std::cout << s.pNucSeq->toString( ) << std::endl; } );

    CppSQLiteExtQueryStatement<int> statement( xTestDB, "SELECT number FROM table1" );
    statement.vExecuteAndForAllRowsUnpackedDo( []( int64_t i ) { std::cout << i << std::endl; } );


} // function

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