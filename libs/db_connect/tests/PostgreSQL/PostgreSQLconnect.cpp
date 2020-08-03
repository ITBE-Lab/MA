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
#include <sql_api.h>
#include <db_con_pool.h>
#include <util/threadPool.h>

// #include "libpq-fe.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef USE_PG
static void exit_nicely( PGconn* conn )
{
    PQfinish( conn );
    exit( 1 );
}

int main_( int argc, char** argv )
{
    const char* conninfo;
    PGconn* conn;
    PGresult* res;
    int nFields;
    int i, j;

    /*
     * If the user supplies a parameter on the command line, use it as the
     * conninfo string; otherwise default to setting dbname=postgres and using
     * environment variables or defaults for all other connection parameters.
     */
    if( argc > 1 )
        conninfo = argv[ 1 ];
    else
        conninfo = "dbname = postgres user = postgres password=admin";

    /* Make a connection to the database */
    conn = PQconnectdb( conninfo );

    /* Check to see that the backend connection was successfully made */
    if( PQstatus( conn ) != CONNECTION_OK )
    {
        fprintf( stderr, "Connection to database failed: %s", PQerrorMessage( conn ) );
        exit_nicely( conn );
    }

    /* Set always-secure search path, so malicious users can't take control. */
    res = PQexec( conn, "SELECT pg_catalog.set_config('search_path', '', false)" );
    if( PQresultStatus( res ) != PGRES_TUPLES_OK )
    {
        fprintf( stderr, "SET failed: %s", PQerrorMessage( conn ) );
        PQclear( res );
        exit_nicely( conn );
    }

    /*
     * Should PQclear PGresult whenever it is no longer needed to avoid memory
     * leaks
     */
    PQclear( res );

    /*
     * Our test case here involves using a cursor, for which we must be inside
     * a transaction block.  We could do the whole thing with a single
     * PQexec() of "select * from pg_database", but that's too trivial to make
     * a good example.
     */

    /* Start a transaction block */
    res = PQexec( conn, "BEGIN" );
    if( PQresultStatus( res ) != PGRES_COMMAND_OK )
    {
        fprintf( stderr, "BEGIN command failed: %s", PQerrorMessage( conn ) );
        PQclear( res );
        exit_nicely( conn );
    }
    PQclear( res );

    /*
     * Fetch rows from pg_database, the system catalog of databases
     */
    res = PQexec( conn, "DECLARE myportal CURSOR FOR select * from pg_database" );
    if( PQresultStatus( res ) != PGRES_COMMAND_OK )
    {
        fprintf( stderr, "DECLARE CURSOR failed: %s", PQerrorMessage( conn ) );
        PQclear( res );
        exit_nicely( conn );
    }
    PQclear( res );

    res = PQexec( conn, "FETCH ALL in myportal" );
    if( PQresultStatus( res ) != PGRES_TUPLES_OK )
    {
        fprintf( stderr, "FETCH ALL failed: %s", PQerrorMessage( conn ) );
        PQclear( res );
        exit_nicely( conn );
    }

    /* first, print out the attribute names */
    nFields = PQnfields( res );
    for( i = 0; i < nFields; i++ )
        printf( "%-15s", PQfname( res, i ) );
    printf( "\n\n" );

    /* next, print out the rows */
    for( i = 0; i < PQntuples( res ); i++ )
    {
        for( j = 0; j < nFields; j++ )
            printf( "%-15s", PQgetvalue( res, i, j ) );
        printf( "\n" );
    }

    PQclear( res );

    /* close the portal ... we don't bother to check for errors ... */
    res = PQexec( conn, "CLOSE myportal" );
    PQclear( res );

    /* end the transaction */
    res = PQexec( conn, "END" );
    PQclear( res );

    /* close the connection to the database and cleanup */
    PQfinish( conn );

    return 0;
}

#if 0
// Example 33.2.libpq Example Program 2


/*
 * src/test/examples/testlibpq2.c
 *
 *
 * testlibpq2.c
 *      Test of the asynchronous notification interface
 *
 * Start this program, then from psql in another window do
 *   NOTIFY TBL2;
 * Repeat four times to get this program to exit.
 *
 * Or, if you want to get fancy, try this:
 * populate a database with the following commands
 * (provided in src/test/examples/testlibpq2.sql):
 *
 *   CREATE SCHEMA TESTLIBPQ2;
 *   SET search_path = TESTLIBPQ2;
 *   CREATE TABLE TBL1 (i int4);
 *   CREATE TABLE TBL2 (i int4);
 *   CREATE RULE r1 AS ON INSERT TO TBL1 DO
 *     (INSERT INTO TBL2 VALUES (new.i); NOTIFY TBL2);
 *
 * Start this program, then from psql do this four times:
 *
 *   INSERT INTO TESTLIBPQ2.TBL1 VALUES (10);
 */

#ifdef WIN32
#include <windows.h>
#endif
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#ifdef HAVE_SYS_SELECT_H
#include <sys/select.h>
#endif

#include "libpq-fe.h"

static void exit_nicely( PGconn* conn )
{
    PQfinish( conn );
    exit( 1 );
}

int main( int argc, char** argv )
{
    const char* conninfo;
    PGconn* conn;
    PGresult* res;
    PGnotify* notify;
    int nnotifies;

    /*
     * If the user supplies a parameter on the command line, use it as the
     * conninfo string; otherwise default to setting dbname=postgres and using
     * environment variables or defaults for all other connection parameters.
     */
    if( argc > 1 )
        conninfo = argv[ 1 ];
    else
        conninfo = "dbname = postgres";

    /* Make a connection to the database */
    conn = PQconnectdb( conninfo );

    /* Check to see that the backend connection was successfully made */
    if( PQstatus( conn ) != CONNECTION_OK )
    {
        fprintf( stderr, "Connection to database failed: %s", PQerrorMessage( conn ) );
        exit_nicely( conn );
    }

    /* Set always-secure search path, so malicious users can't take control. */
    res = PQexec( conn, "SELECT pg_catalog.set_config('search_path', '', false)" );
    if( PQresultStatus( res ) != PGRES_TUPLES_OK )
    {
        fprintf( stderr, "SET failed: %s", PQerrorMessage( conn ) );
        PQclear( res );
        exit_nicely( conn );
    }

    /*
     * Should PQclear PGresult whenever it is no longer needed to avoid memory
     * leaks
     */
    PQclear( res );

    /*
     * Issue LISTEN command to enable notifications from the rule's NOTIFY.
     */
    res = PQexec( conn, "LISTEN TBL2" );
    if( PQresultStatus( res ) != PGRES_COMMAND_OK )
    {
        fprintf( stderr, "LISTEN command failed: %s", PQerrorMessage( conn ) );
        PQclear( res );
        exit_nicely( conn );
    }
    PQclear( res );

    /* Quit after four notifies are received. */
    nnotifies = 0;
    while( nnotifies < 4 )
    {
        /*
         * Sleep until something happens on the connection.  We use select(2)
         * to wait for input, but you could also use poll() or similar
         * facilities.
         */
        int sock;
        fd_set input_mask;

        sock = PQsocket( conn );

        if( sock < 0 )
            break; /* shouldn't happen */

        FD_ZERO( &input_mask );
        FD_SET( sock, &input_mask );

        if( select( sock + 1, &input_mask, NULL, NULL, NULL ) < 0 )
        {
            fprintf( stderr, "select() failed: %s\n", strerror( errno ) );
            exit_nicely( conn );
        }

        /* Now check for input */
        PQconsumeInput( conn );
        while( ( notify = PQnotifies( conn ) ) != NULL )
        {
            fprintf( stderr, "ASYNC NOTIFY of '%s' received from backend PID %d\n", notify->relname, notify->be_pid );
            PQfreemem( notify );
            nnotifies++;
            PQconsumeInput( conn );
        }
    }

    fprintf( stderr, "Done.\n" );

    /* close the connection to the database and cleanup */
    PQfinish( conn );

    return 0;
}


Example 33.3.libpq Example Program 3


/*
 * src/test/examples/testlibpq3.c
 *
 *
 * testlibpq3.c
 *      Test out-of-line parameters and binary I/O.
 *
 * Before running this, populate a database with the following commands
 * (provided in src/test/examples/testlibpq3.sql):
 *
 * CREATE SCHEMA testlibpq3;
 * SET search_path = testlibpq3;
 * CREATE TABLE test1 (i int4, t text, b bytea);
 * INSERT INTO test1 values (1, 'joe''s place', '\\000\\001\\002\\003\\004');
 * INSERT INTO test1 values (2, 'ho there', '\\004\\003\\002\\001\\000');
 *
 * The expected output is:
 *
 * tuple 0: got
 *  i = (4 bytes) 1
 *  t = (11 bytes) 'joe's place'
 *  b = (5 bytes) \000\001\002\003\004
 *
 * tuple 0: got
 *  i = (4 bytes) 2
 *  t = (8 bytes) 'ho there'
 *  b = (5 bytes) \004\003\002\001\000
 */

#ifdef WIN32
#include <windows.h>
#endif

#include "libpq-fe.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

/* for ntohl/htonl */
#include <arpa/inet.h>
#include <netinet/in.h>


    static void
    exit_nicely( PGconn* conn )
{
    PQfinish( conn );
    exit( 1 );
}

/*
 * This function prints a query result that is a binary-format fetch from
 * a table defined as in the comment above.  We split it out because the
 * main() function uses it twice.
 */
static void show_binary_results( PGresult* res )
{
    int i, j;
    int i_fnum, t_fnum, b_fnum;

    /* Use PQfnumber to avoid assumptions about field order in result */
    i_fnum = PQfnumber( res, "i" );
    t_fnum = PQfnumber( res, "t" );
    b_fnum = PQfnumber( res, "b" );

    for( i = 0; i < PQntuples( res ); i++ )
    {
        char* iptr;
        char* tptr;
        char* bptr;
        int blen;
        int ival;

        /* Get the field values (we ignore possibility they are null!) */
        iptr = PQgetvalue( res, i, i_fnum );
        tptr = PQgetvalue( res, i, t_fnum );
        bptr = PQgetvalue( res, i, b_fnum );

        /*
         * The binary representation of INT4 is in network byte order, which
         * we'd better coerce to the local byte order.
         */
        ival = ntohl( *( (uint32_t*)iptr ) );

        /*
         * The binary representation of TEXT is, well, text, and since libpq
         * was nice enough to append a zero byte to it, it'll work just fine
         * as a C string.
         *
         * The binary representation of BYTEA is a bunch of bytes, which could
         * include embedded nulls so we have to pay attention to field length.
         */
        blen = PQgetlength( res, i, b_fnum );

        printf( "tuple %d: got\n", i );
        printf( " i = (%d bytes) %d\n", PQgetlength( res, i, i_fnum ), ival );
        printf( " t = (%d bytes) '%s'\n", PQgetlength( res, i, t_fnum ), tptr );
        printf( " b = (%d bytes) ", blen );
        for( j = 0; j < blen; j++ )
            printf( "\\%03o", bptr[ j ] );
        printf( "\n\n" );
    }
}

int main( int argc, char** argv )
{
    const char* conninfo;
    PGconn* conn;
    PGresult* res;
    const char* paramValues[ 1 ];
    int paramLengths[ 1 ];
    int paramFormats[ 1 ];
    uint32_t binaryIntVal;

    /*
     * If the user supplies a parameter on the command line, use it as the
     * conninfo string; otherwise default to setting dbname=postgres and using
     * environment variables or defaults for all other connection parameters.
     */
    if( argc > 1 )
        conninfo = argv[ 1 ];
    else
        conninfo = "dbname = postgres";

    /* Make a connection to the database */
    conn = PQconnectdb( conninfo );

    /* Check to see that the backend connection was successfully made */
    if( PQstatus( conn ) != CONNECTION_OK )
    {
        fprintf( stderr, "Connection to database failed: %s", PQerrorMessage( conn ) );
        exit_nicely( conn );
    }

    /* Set always-secure search path, so malicious users can't take control. */
    res = PQexec( conn, "SET search_path = testlibpq3" );
    if( PQresultStatus( res ) != PGRES_COMMAND_OK )
    {
        fprintf( stderr, "SET failed: %s", PQerrorMessage( conn ) );
        PQclear( res );
        exit_nicely( conn );
    }
    PQclear( res );

    /*
     * The point of this program is to illustrate use of PQexecParams() with
     * out-of-line parameters, as well as binary transmission of data.
     *
     * This first example transmits the parameters as text, but receives the
     * results in binary format.  By using out-of-line parameters we can avoid
     * a lot of tedious mucking about with quoting and escaping, even though
     * the data is text.  Notice how we don't have to do anything special with
     * the quote mark in the parameter value.
     */

    /* Here is our out-of-line parameter value */
    paramValues[ 0 ] = "joe's place";

    res = PQexecParams( conn,
                        "SELECT * FROM test1 WHERE t = $1",
                        1, /* one param */
                        NULL, /* let the backend deduce param type */
                        paramValues,
                        NULL, /* don't need param lengths since text */
                        NULL, /* default to all text params */
                        1 ); /* ask for binary results */

    if( PQresultStatus( res ) != PGRES_TUPLES_OK )
    {
        fprintf( stderr, "SELECT failed: %s", PQerrorMessage( conn ) );
        PQclear( res );
        exit_nicely( conn );
    }

    show_binary_results( res );

    PQclear( res );

    /*
     * In this second example we transmit an integer parameter in binary form,
     * and again retrieve the results in binary form.
     *
     * Although we tell PQexecParams we are letting the backend deduce
     * parameter type, we really force the decision by casting the parameter
     * symbol in the query text.  This is a good safety measure when sending
     * binary parameters.
     */

    /* Convert integer value "2" to network byte order */
    binaryIntVal = htonl( (uint32_t)2 );

    /* Set up parameter arrays for PQexecParams */
    paramValues[ 0 ] = (char*)&binaryIntVal;
    paramLengths[ 0 ] = sizeof( binaryIntVal );
    paramFormats[ 0 ] = 1; /* binary */

    res = PQexecParams( conn,
                        "SELECT * FROM test1 WHERE i = $1::int4",
                        1, /* one param */
                        NULL, /* let the backend deduce param type */
                        paramValues,
                        paramLengths,
                        paramFormats,
                        1 ); /* ask for binary results */

    if( PQresultStatus( res ) != PGRES_TUPLES_OK )
    {
        fprintf( stderr, "SELECT failed: %s", PQerrorMessage( conn ) );
        PQclear( res );
        exit_nicely( conn );
    }

    show_binary_results( res );

    PQclear( res );

    /* close the connection to the database and cleanup */
    PQfinish( conn );

    return 0;
}
#endif
#endif

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


int main( int argc, char** argv )
{
    try
    {
        // definition of database connection in json:
        auto jDBConfig = json{
            { SCHEMA, { { NAME, "sv_db" } /* , { FLAGS, { DROP_ON_CLOSURE } } */ } },
        // { TEMPORARY, true },
#ifdef USE_PG
            { CONNECTION, { { HOSTNAME, "localhost" }, { USER, "postgres" }, { PASSWORD, "admin" }, { PORT, 0 } } }
#endif
#ifdef USE_MSQL
            { CONNECTION, { { HOSTNAME, "localhost" }, { USER, "root" }, { PASSWORD, "admin" }, { PORT, 0 } } }
#endif
        };

        std::vector<std::future<void>> vFutures;
        {
            SQLDBConPool<DBConn> xDBPool( 1, jDBConfig );
            for( int i = 0; i < 1 /* 20 */; i++ )
                // type behind auto: std::shared_ptr<SQLDBConPool<MySQLConDB>::PooledSQLDBCon>
                vFutures.push_back( xDBPool.enqueue( []( auto pDBCon ) {
                    doNoExcept(
                        [ & ] {
                            // typedef decltype( *pDBCon ) Type;
                            // typedef typename SubType :: element_type Type;
                            // using Type = typename decltype( pDBCon )::element_type;
                            std::cout << "Job executed in task: " << pDBCon->getTaskId( ) << std::endl;

                            // Database Test
                            checkDB( pDBCon, pDBCon->getTaskId( ) );

                            pDBCon->doPoolSafe( [] { std::cout << "This print is pool safe ..." << std::endl; } );
                        },
                        "Problem during thread execution" );
                } ) );
        } // close the pool

        // Get all future exception safe
        for( auto& rFurture : vFutures )
            doNoExcept( [ & ] { rFurture.get( ); } );

        std::cout << "ALL WORK DONE ..." << std::endl;
        // #ifdef _MSC_VER
        //         int i;
        //         std::cin >> i;
        // #endif
        return 0;

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
} // main
