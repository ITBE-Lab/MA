/* Authors: Arne Kutzner and Markus Schmidt
 * Created: Jan. 2020
 * MIT License
 * @file common.h
 * @brief General API for SQL-like databases. Use this API in combination with a connector to work with SQL databases.
 */

#pragma once
// #define SQL_VERBOSE // define me if required

#ifdef _MSC_VER
#pragma warning( disable : 4996 ) // suppress warnings regarding the use of strerror
#endif

#include <cerrno> // error management for file I/O
#include <cstring> // error management for file I/O
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <vector>

// For getting the code working with gcc 6.x.x compiler
#if( __GNUC__ && ( __GNUC__ < 7 ) )
#include <experimental/tuple>
#define STD_APPLY std::experimental::apply
#else
#include <tuple>
#define STD_APPLY std::apply
#endif

#if defined( __GNUC__ ) && ( __GNUC__ < 8 )
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

// Import JSON support
#include <json.hpp>
using nlohmann::json;

/** @brief Returns a string that contains info regarding types and values of an argument pack.
 *  Primarily for debugging purposes...
 */
template <typename... ArgTypes> std::string dumpArgPack( ArgTypes&&... args )
{
    std::stringstream xInfStream;
    ( ( xInfStream << typeid( args ).name( ) << ":" << args << "  " ), ... );
    return xInfStream.str( );
} // helper function

/** @brief Delivers hex string for integer value.
 *  See: https://stackoverflow.com/questions/5100718/integer-to-hex-string-in-c
 */
template <typename T> inline std::string intToHex( T val, size_t width = sizeof( T ) * 2 )
{
    std::stringstream ss;
    ss << std::setfill( '0' ) << std::setw( width ) << std::hex << ( val | 0 );
    return ss.str( );
} // helper function

/** @brief Executes the func so that exceptions are swallowed.
 *         Returns true, if an exceptions was dropped.
 *         Returns false, if execution was OK.
 *  @detail For destructors threads and tests so that they do not throw.
 */
template <typename F> inline bool doNoExcept( F&& func, std::string sInfo = "" )
{
    try
    {
        func( );
        return false;
    } // try
    catch( std::exception& rxExcept )
    {
        std::cout << sInfo << ( sInfo.empty( ) ? "" : " - " )
                  << std::string( "Dropped exception:\n" ) + rxExcept.what( ) << std::endl;
        return true;
    } // catch
    catch( ... )
    {
        std::cout << sInfo << ( sInfo.empty( ) ? "" : "\n" ) << "Dropped unknown exception." << std::endl;
        return true;
    } // catch
} // method


/** @brief An instance represents a connection to DB-system like MySQL etc.
 *   DBImpl must be a class that implements a database API.
 *  @detail A single connection is supposed not to be thread-safe. For multiple threads use a connection pool.
 */
template <typename DBImpl> class SQLDB : public DBImpl
{
    /** brief The guard guarantees that the transaction commits of the object runs out of scope.
     */
    class GuardedTransaction
    {
        SQLDB<DBImpl>& rHost; // host DB connection
        std::shared_ptr<bool> pHostTombStone; // for checking the aliveness of the host
        bool bStateCommitted; // Is true if the transaction has been committed, false otherwise

      public:
        GuardedTransaction( const GuardedTransaction& ) = delete; // no transaction copies

        GuardedTransaction( SQLDB<DBImpl>& rHost, std::shared_ptr<bool> pHostTombStone )
            : rHost( rHost ), pHostTombStone( pHostTombStone ), bStateCommitted( true )
        {
            this->start( );
        } // constructor

        /** brief Commit the transaction */
        void commit( )
        {
            if( bStateCommitted )
                throw std::runtime_error( "SQL-DB Transaction:\nYou tried to commit after committing already." );
            if( *pHostTombStone )
                throw std::runtime_error( "SQL-DB Transaction:\nYou tried to commit although the database connection "
                                          "has been gone already." );
            rHost.commitTrxn( );
            this->bStateCommitted = true;
        } // method

        /** brief Start the transaction */
        void start( )
        {
            if( !bStateCommitted )
                throw std::runtime_error( "SQL-DB Transaction:\nYou have to commit before starting a transaction." );
            if( *pHostTombStone )
                throw std::runtime_error(
                    "SQL-DB Transaction:\nYou tried to start a transaction although the database connection "
                    "has been gone already." );
            rHost.startTrxn( );
            this->bStateCommitted = false;
        } // method

        /* Commit the transaction if the corresponding about runs out of scope.*/
        ~GuardedTransaction( )
        {
            this->commit( );
        } // destructor
    }; // class (GuardedTransaction)

    std::shared_ptr<bool> pTombStone; // So long the database connection is alive the stone says false.
                                      // As soon as the DB is gone, the stone says true.


  public:
    typedef std::unique_ptr<GuardedTransaction> uniqueGuardedTrxnType;
    typedef std::shared_ptr<GuardedTransaction> sharedGuardedTrxnType;
    typedef DBImpl DBImplForward; // Forwarding of the template parameter type
    const std::string sConId; // unique connection ID with respect to the current machine
    const std::string sSchemaName; // unique connection ID with respect to the current machine
    const bool bTemporary; // this is true the schema shall self delete in its destructor

    SQLDB( const SQLDB& ) = delete; // no DB connection copies

    // For serializing all table insertions with respect to the current DB-connection.
    // Necessary for getting the correct last insertion ID in case of an auto-increment column. (e.g. with MySQL)
    // FIXME: Remove the lock, because it is not necessary.
    std::mutex pGlobalInsertLock;

    SQLDB( const json& jDBConData = json{} )
        : DBImpl( jDBConData ),
          pTombStone( std::make_shared<bool>( false ) ), // initialize tombstone
          sConId( intToHex( reinterpret_cast<uint64_t>( this ) ) ), // use the address for id creation
          sSchemaName( jDBConData.count( "SCHEMA" ) ? jDBConData[ "SCHEMA" ] : "" ),
          bTemporary( jDBConData.count( "TEMPORARY" ) == 1 && jDBConData[ "TEMPORARY" ].get<bool>() )
    {} // constructor

    /**
     * @brief initialize DB connection from schema name
     * @details
     * this constructor is exported to python
     */
    SQLDB( std::string sSchemaName ) : SQLDB( json{{"SCHEMA", sSchemaName}} )
    {}

    ~SQLDB( )
    {
        if( bTemporary && !sSchemaName.empty( ) )
            dropSchema( sSchemaName );
        // Inform about the "death" of the database connection.
        *pTombStone = true;
        // Now the destructor of the DB implementor is called ...
    }; // destructor

    int64_t lastRowId( )
    {
        throw std::runtime_error( "lastRowId of SQLDB is not implemented yet." );
    } // method

    /** @brief Directly execute the statement in DB. */
    void execStmt( const std::string& rsStmtTxt )
    {
        DBImpl::execSQL( rsStmtTxt );
    } // method

    /** @brief A commitGuard immediately starts a transaction and commits it, if it goes out of scope.
     * Destructing the guard after destructing the database leads to a crash.
     */
    std::unique_ptr<GuardedTransaction> uniqueGuardedTrxn( )
    {
        return std::make_unique<GuardedTransaction>( *this, pTombStone );
    } // method

    std::shared_ptr<GuardedTransaction> sharedGuardedTrxn( )
    {
        return std::make_shared<GuardedTransaction>( *this, pTombStone );
    } // method

    // WARNING: If indexExists is coupled with an index creation, it should be executed "pool safe".
    bool indexExists( const std::string& rsTblName, const std::string& rsIdxName )
    {
        return DBImpl::indexExistsInDB( rsTblName, rsIdxName );
    } // method

    /** @brief Drop Database Schema */
    void dropSchema( const std::string& rsSchemaName )
    {
        DBImpl::execSQL( "DROP DATABASE " + rsSchemaName );
    } // method

    /** @brief Use Database Schema */
    void useSchema( const std::string& rsSchemaName )
    {
        DBImpl::execSQL( "USE DATABASE " + rsSchemaName );
    } // method

    /** @brief This function should be redefined in pooled SQL connections so that the function is executed in mutex
     * protected environment. In a single connection we simply execute the function.
     */
    template <typename F> void doPoolSafe( F&& func )
    {
        func( );
    } // method
}; // SQLDB


template <typename DBCon> class SQLStatement
{
  protected:
    std::shared_ptr<DBCon> pDB; // Database connection //--
    std::unique_ptr<typename DBCon::PreparedStmt> pStmt; // pointer to insert statement

  public:
    SQLStatement( const SQLStatement& ) = delete; // no statement copies

    SQLStatement( std::shared_ptr<DBCon> pDB, const std::string& rsStmtText ) //--
        : pDB( pDB ), pStmt( std::make_unique<typename DBCon::PreparedStmt>( this->pDB, rsStmtText ) )
    {} // constructor

    template <typename... ArgTypes> inline int exec( ArgTypes&&... args )
    {
        return static_cast<int>( pStmt->bindAndExec( std::forward<ArgTypes>( args )... ) );
    } // method
}; // class


/* @brief An instance of this class models an SQL query. */
template <typename DBCon, typename... ColTypes> class SQLQuery
{
  private:
    std::shared_ptr<DBCon> pDB; // Database connection //--
    std::unique_ptr<typename DBCon::template PreparedQuery<ColTypes...>> pQuery; // pointer to insert statement
    const std::string sStmtText; // backup of the query statement text. (used for verbose error reporting)
    /* Deprecated. Execute and bind the query args., but do not fetch the first row. */

    template <typename F> void forAllRowsDo( F&& func )
    {
        // Fetch all rows and apply func to all of them
        while( pQuery->fetchNextRow( ) )
            // unpack the tuple and call func (C++17 fold expression)
            STD_APPLY( func, pQuery->tCellValues );
    } // method

  public:
    /** @brief Constructs a query instance using the statement rsStmtText */
    SQLQuery( std::shared_ptr<DBCon> pDB, const std::string& rsStmtText )
        : pDB( pDB ),
          pQuery( std::make_unique<typename DBCon::template PreparedQuery<ColTypes...>>( this->pDB, rsStmtText ) ),
          sStmtText( rsStmtText )
    {} // constructor

    /** @brief Execute the query */
    template <typename... ArgTypes> void exec( ArgTypes&&... args )
    {
        // Execute statement
        pQuery->bindAndExec( std::forward<ArgTypes>( args )... );
        // Bind outcome of statement execution for later fetching
        pQuery->bindResult( );
    } // method

    /** @brief Execute the query and fetch first row */
    template <typename... ArgTypes> bool execAndFetch( ArgTypes&&... args )
    {
        // Execute the query
        this->exec( std::forward<ArgTypes>( args )... );
        // Fetch first row for getting correct EOF value.
        return this->next( );
    } // method

    /** @brief Execute the query passing args and call the function func for all rows of the query outcome. */
    template <typename F, typename... ArgTypes> void execAndForAll( F&& func, ArgTypes&&... args )
    {
        this->exec( std::forward<ArgTypes>( args )... );
        this->forAllRowsDo( std::forward<F>( func ) );
    } // method

    /** @brief Execute the query and get the n-th cell of the first row of a query. */
    template <int COL_NUM, typename... ArgTypes>
    typename std::tuple_element<COL_NUM, std::tuple<ColTypes...>>::type execAndGetNthCell( ArgTypes&&... args )
    {
        // Execute query, bind, and fetch first row
        this->exec( std::forward<ArgTypes>( args )... );
        if( !pQuery->fetchNextRow( ) )
            throw std::runtime_error( "SQLQuery: Tried to retrieve first row for an empty query outcome. For stmt:" +
                                      sStmtText + "\nArgs: " + dumpArgPack( args... ) );
        return std::get<COL_NUM>( pQuery->tCellValues );
    } // method

    /** @brief Executes the query and delivers the first value of the first row as result.
     *  WARNING: The type is automatically the type of the first column. (FIXME: fix this!)
     */
    typedef typename std::tuple_element<0, std::tuple<ColTypes...>>::type FstColType;
    template <typename... ArgTypes> FstColType scalar( ArgTypes&&... args )
    {
        return execAndGetNthCell<0>( std::forward<ArgTypes>( args )... );
    } // method

    /** @brief Delivers a reference to a tuple holding the data of the current row.
     * WARNING: For efficiency reasons, the tuple is delivered by reference. The returned reference
     * stays valid for the lifetime of the current query object.
     * https://stackoverflow.com/questions/24725740/using-auto-and-decltype-to-return-reference-from-function-in-templated-class
     */
    inline auto get( ) -> decltype( pQuery->getCellValues( ) )&
    {
        try
        {
            return pQuery->getCellValues( );
        } // try
        catch( const std::exception& rxExept )
        {
            throw std::runtime_error( "SQLQuery: next() received an internal error for stmt: " + sStmtText +
                                      "\nInternal Error: " + rxExept.what( ) );
        } // catch
    } // method

    /** @brief Fetch next row of query.
     *  Returns true, if a row could be fetched successfully, false otherwise (indicating EOF).
     *	The fresh row can be received via get.
     */
    inline bool next( )
    {
        return pQuery->fetchNextRow( );
    } // method

    /** @brief End Of File. Delivers true, if the previous call of next() did not deliver a fresh row.
     *  This kind of eof works in a "look back" fashion; it tells if the previous fetch failed or not.
     */
    inline bool eof( )
    {
        try
        {
            return pQuery->eofLookBack( );
        } // try
        catch( const std::exception& rxExept )
        {
            throw std::runtime_error( "SQLQuery: eof() received an internal error for stmt: " + sStmtText +
                                      "\nInternal Error: " + rxExept.what( ) );
        } // catch
    } // method

    /* Delivers the outcome of a query as list of tuples.
     * DEPRECATED!
     */
    template <class... ArgTypes>
    std::vector<typename std::tuple<ColTypes...>> executeAndStoreAllInVector( ArgTypes&&... args )
    {
        // throw std::runtime_error( "executeAndStoreAllInVector is untested." );
        // Returned vector that receives column of table
        std::vector<typename std::tuple<ColTypes...>> vRetVec;
        // Execute query, bind, and fetch first row
        this->exec( std::forward<ArgTypes>( args )... );
        // Fetch all rows and add them to the returned vector
        while( pQuery->fetchNextRow( ) )
            vRetVec.push_back( pQuery->tCellValues );
        return vRetVec;
    } // template method


    /* Delivers one column of the outcome of query as vector.
     * DEPRECATED!
     */
    template <size_t COL_NUM, class... ArgTypes>
    std::vector<typename std::tuple_element<COL_NUM, std::tuple<ColTypes...>>::type>
    executeAndStoreInVector( ArgTypes&&... args )
    {
        // throw std::runtime_error( "executeAndStoreAllInVector is untested." );
        // Returned vector that receives column of table
        std::vector<typename std::tuple_element<COL_NUM, std::tuple<ColTypes...>>::type> vRetVec;
        // Execute query, bind, and fetch first row
        this->exec( std::forward<ArgTypes>( args )... );
        // Fetch all rows and add them to the returned vector
        while( pQuery->fetchNextRow( ) )
            vRetVec.emplace_back( std::get<COL_NUM>( pQuery->tCellValues ) );
        return vRetVec;
    } // template method
}; // class


// Constants for table definitions via json
const std::string TABLE_NAME = "TABLE_NAME";
const std::string TABLE_COLUMNS = "TABLE_COLUMNS";
const std::string SQLITE_EXTRA = "SQLITE_EXTRA";
const std::string COLUMN_NAME = "COLUMN_NAME";
const std::string CONSTRAINTS = "CONSTRAINTS";
const std::string PLACEHOLDER = "PLACEHOLDER";
const std::string CPP_EXTRA = "CPP_EXTRA";
const std::string SQL_EXTRA = "SQL_EXTRA";
const std::string PRIMARY_KEY = "PRIMARY_KEY";
const std::string FOREIGN_KEY = "FOREIGN_KEY";
const std::string REFERENCES = "REFERENCES";

// Constants for index definitions via json
const std::string INDEX_NAME = "INDEX_NAME";
const std::string INDEX_COLUMNS = "INDEX_COLUMNS";
const std::string WHERE = "WHERE";

// generated columns; see: https://dev.mysql.com/doc/refman/5.7/en/create-table-generated-columns.html
// @todo talk with arne about design here
// part of ColTypes?:
// yes: problem with inserters
// no: TYPE keyword required... (currently this is implemented)
const std::string GENERATED_COLUMNS = "GENERATED_COLUMNS";
// the following keywords do only have an effect on generated columns
const std::string TYPE = "TYPE";
const std::string AS = "AS";
const std::string STORED = "STORED";


/** @brief Serializes a value for CSV representation.
 *  The translation is required by the SQLFileBulkInserter.
 *  NULL values are represented by '\N'.
 */
template <typename DBCon, typename Type> inline std::string csvPrint( DBCon& rxDBCon, Type& rVal )
{
    return "'" + std::to_string( rVal ) + "'";
}; // method

template <typename DBCon> inline std::string csvPrint( DBCon& rxDBCon, const std::string& rVal )
{
    return "'" + rVal + "'";
}; // method

template <typename DBCon> inline std::string csvPrint( DBCon& rxDBCon, const std::nullptr_t& rVal )
{
    return "\\N";
}; // method

/**
 * @brief completely empty struct
 * @details
 * can be used to hold and then extract a parameter pack from a variable.
 * @todo pack -> Pack
 */
template <typename... Args> struct pack
{};

/* Tips: https://stackoverflow.com/questions/667643/mysql-batch-updates-in-c
 * Change the name to SQLTableView
 */
template <typename DBCon, typename... ColTypes> class SQLTable
{
  protected:
    /** @brief: Inner class for representing view to an index.
     * { "INDEX_NAME", "name_of_your_index" } - item is option
     * { "INDEX_COLUMNS", "columns in comma separated form "} - duty
     * { "WHERE", "predicate expression" } - item is option
     */
    class SQLIndexView
    {
        const json jIndexDef; // copy of JSON index definition

        std::string makeIndexCreateStmt( const std::string& rsTblName, const std::string& rsIdxName )
        {
            if( jIndexDef.count( INDEX_COLUMNS ) != 1 )
                throw std::runtime_error(
                    "The JSON definition for a SQL index requires one INDEX_COLUMNS item exactly." );
            std::string sCols = jIndexDef[ INDEX_COLUMNS ];
            std::string sStmt = "CREATE INDEX " + rsIdxName + " ON " + rsTblName + "(" + sCols + ")";

            // WHERE is not supported by MySQL
            if( jIndexDef.count( WHERE ) )
            {
                if( DBCon::supportsWhereClauseForIndex( ) )
                    sStmt.append( " WHERE " ).append( jIndexDef[ WHERE ] );
                else
                    std::cout << "WARNING: WHERE clause ignored." << std::endl;
            } // if
#ifdef SQL_VERBOSE
            std::cout << "Index creation statement: " << sStmt << std::endl;
#endif
            return sStmt;
        } // method

      public:
        SQLIndexView( const SQLTable<DBCon, ColTypes...>* pTable,
                      const json& rjIndexDef ) // index definition in JSON)
            : jIndexDef( rjIndexDef ) // copy the index def. to attribute
        {
            // If there is no index name in the JSON, generate one.
            std::string sIdxName = jIndexDef.count( INDEX_NAME ) ? std::string( jIndexDef[ INDEX_NAME ] )
                                                                 : pTable->getTableName( ) + "_idx";

            // In a pooled environment the creation of indexes should be serialzed.
            pTable->pDB->doPoolSafe( [&] {
                // When we check the existence of the index as well as during creation,
                // we require an exclusive lock on the database connection.
                if( !pTable->pDB->indexExists( pTable->getTableName( ), sIdxName ) )
                {
                    pTable->pDB->execSQL( makeIndexCreateStmt( pTable->getTableName( ), sIdxName ) );
                } // if
                else
                {
#ifdef SQL_VERBOSE
                    std::cout << "Index exists already, skip creation." << std::endl;
#endif
                }
            } ); // doPoolSafe
        } // constructor
    }; // inner class SQL Index

  public:
    /** @brief: Implements the concept of bulk inserts for table views.
     *  The bulk inserter always inserts NULL's on columns having type std::nullptr_t.
     */
    template <size_t BUF_SIZE,
              typename... InsTypes> // InsTypes - either equal to corresponding type in ColTypes or std::nullptr_t
    class SQLBulkInserter
    {
      private:
        std::array<std::tuple<InsTypes...>, BUF_SIZE> aBuf; // buffer of values inserted by a single op.
        size_t uiInsPos; // index of the next insert into the array
        std::unique_ptr<typename DBCon::PreparedStmt> pBulkInsertStmt; // pointer to prepared bulk insert statement
        std::unique_ptr<typename DBCon::PreparedStmt> pSingleInsertStmt; // required for buffer flush

#if BULK_INSERT_BY_TUPLE_CATENATION // Nice Haskell like approach, but apparently to much for most C++ compilers ...
        // Implementation part of cat_arr_of_tpl
        template <typename Array, std::size_t... I>
        auto cat_arr_of_tpl_impl( const Array& a, std::index_sequence<I...> )
        {
            return std::tuple_cat( std::move( a[ I ] )... );
        } // meta

        /** @brief Concatenates an array of tuples and returns the resulting single tuple.
         * IMPORTANT:
         * Tuple elements are moved for efficiency, so the array and its elements should not be moved after
         * catenation.
         */
        template <typename T, std::size_t N, typename Indices = std::make_index_sequence<N>>
        auto cat_arr_of_tpl( const std::array<T, N>& a )
        {
            return cat_arr_of_tpl_impl( a, Indices{} );
        } // meta

        /** @brief Inserts the buffer content into the table by executing a bulk insert statement.
         * Design alternative: iterate over the array and bind each tuple separately.
         * Discussion: The latter approach requires more code but avoids template instantiation problems.
         */
        inline void bulkInsert( )
        {
            auto tCatTpl = cat_arr_of_tpl( aBuf ); // all tuples of aBuf concatenated to a single one
            // DEBUG: std::cout << typeid(tCatTpl).name() << std::endl;

            // std::apply transforms the concatenated tuple to a argument pack and passes this to the lambda.
            // The lambda forwards the argument to bindAndExec via references for avoiding copies.
            // This statement creates trouble for some GCC compilers (template instantiation depth exceeds maximum)
            // ...
            STD_APPLY( [&]( auto&... args ) { pBulkInsertStmt->bindAndExec( args... ); }, tCatTpl );
        } // method
#else // Iterating approach for parameter binding with bulk inserts. (Lacks beauty, but better manageable for compilers)

        /** @brief Unpack the tuple to an argument pack and bind using this pack.
         *  (Part of bulkInsert)
         */
        template <int OFFSET, typename... Args> void doSingleBind( const std::tuple<Args...>& rTpl )
        {
            STD_APPLY( [&]( auto&... args ) //
                       { pBulkInsertStmt->template bind<OFFSET>( args... ); },
                       rTpl );
        } // method

        /** @brief Implementation part of forAllDoBind */
        template <typename T, std::size_t N, std::size_t... Idx>
        void forAllDoBindImpl( const std::array<T, N>& a, std::index_sequence<Idx...> )
        {
            // Requires C++17 folded expression support ...
            ( doSingleBind<Idx>( a[ Idx ] ), ... );
        } // meta

        /** @brief Compile time iteration over the array a for binding each tuple in the array (each row) */
        template <typename T, std::size_t N> void forAllDoBind( const std::array<T, N>& a )
        {
            forAllDoBindImpl( a, std::make_index_sequence<N>{} );
        } // meta

        /** @brief Inserts the buffer content into the table by executing a bulk insert statement.
         *  Design alternative: iterate over the array and bind each tuple separately.
         *  Discussion: The latter approach requires more code but avoids the template instantiation problems.
         */
        inline void bulkInsert( )
        {
            // Iterates over the array during compile time.
            forAllDoBind( aBuf );
            // After binding all arguments, actually execute the insert statement.
            pBulkInsertStmt->exec( );
        } // method
#endif // alternative approach

      public:
        /** @brief Construct bulk-inserter.
         *  The table must be part of the database already, or the construction fails.
         *  Implementation detail:
         *  Because we can not guarantee the life of host table for all of the life of the bulk inserter,
         *  we should not keep a reference to the table. Instead we copy the required info of the table.
         */
        SQLBulkInserter( const SQLTable& rxHost )
            : uiInsPos( 0 ),
              pBulkInsertStmt( // compile bulk insert stmt
                  std::make_unique<typename DBCon::PreparedStmt>( rxHost.pDB, rxHost.makeInsertStmt( BUF_SIZE ) ) ),
              pSingleInsertStmt( // compile single insert stmt
                  std::make_unique<typename DBCon::PreparedStmt>( rxHost.pDB, rxHost.makeInsertStmt( 1 ) ) )
        {
            static_assert( sizeof...( InsTypes ) == sizeof...( ColTypes ) );
        } // constructor

        /** @brief Inserts a row into the table via a bulk-insert approach.
         * Reasonable addition: insert using moves.
         */
        inline void insert( const InsTypes&... args )
        {
            aBuf[ uiInsPos ] = std::tuple<InsTypes...>( args... ); // This triggers a copy...
            if( ++uiInsPos == BUF_SIZE )
            {
                // bulkInsert( ) does not inspects uiInsPos anymore. However, it can throw an exceptions.
                // In order to avoid the crash scenario described in flush(),  uiInsPos must be set zero before
                // calling bulkInsert().
                uiInsPos = 0; // reset counter
                bulkInsert( ); // write full buffer to DB
            } // if
        } // method

        /** @brief Flush a non-full buffer to the database.
         */
        inline void flush( )
        {
            // uiInsPos must be set zero before calling bindAndExec, because bindAndExec could throw an exception
            // which in turn triggers a destructor call. However, in the case of an exception, the flush in the
            // destructor should not lead to an additional call of bindAndExec, because it would lead to a crash.
            auto uiInsPosBackup = uiInsPos;
            uiInsPos = 0; // reset insert position counter
            for( size_t uiCount = 0; uiCount < uiInsPosBackup; uiCount++ )
                // Write the current tuple in the array to the DB
                STD_APPLY( [&]( auto&... args ) { pSingleInsertStmt->bindAndExec( args... ); }, aBuf[ uiCount ] );
        } // method

        /** @brief Destructor flushes the buffer.
         */
        ~SQLBulkInserter( )
        {
            // Throwing an exception in a destructor results in undefined behavior.
            // Therefore, we swallow these exceptions and report them via std:cerr.
            doNoExcept( [this] { this->flush( ); } );
        } // destructor
    }; // class


    /** @brief Like a standard SQLBulkInserter, but for the streaming a CSV-file is used in between. This approach
     * is sometimes faster than the standard SQLBulkInserter. IMPORTANT NOTICE: Call release() before a
     * SQLFileBulkInserter runs out of scope. In this case you get an exception, if something goes wrong. If you let
     * do the release-job automatically by the destructor all exception are swallowed and you won't get any feedback
     * about problems.
     */
    template <size_t BUF_SIZE,
              typename... InsTypes> // InsTypes - either equal to corresponding type in ColTypes or std::nullptr_t
    class SQLFileBulkInserter
    {
      private:
        std::array<std::tuple<InsTypes...>, BUF_SIZE> aBuf; // buffer of values inserted by a single op.
        size_t uiInsPos; // index of the next insert into the array
        size_t uiInsCnt; // counter for all inserts so far

        // Tuple printer for ostream.
        // Found at: https://en.cppreference.com/w/cpp/utility/integer_sequence
        template <class Ch, class Tr, class Tuple, std::size_t... Is>
        void prtTpltoStreamImpl( std::basic_ostream<Ch, Tr>& os, const Tuple& t, std::index_sequence<Is...> )
        {
            // Requires C++17 folded expression support ...
            ( ( os << ( Is == 0 ? "" : "\t" ) << csvPrint( *pHostTblDB, std::get<Is>( t ) ) ), ... );
        } // meta

        template <class Ch, class Tr, class... Args>
        void prtTpltoStream( std::basic_ostream<Ch, Tr>& os, const std::tuple<Args...>& t )
        {
            prtTpltoStreamImpl( os, t, std::index_sequence_for<Args...>{} );
            os << "\n";
        } // meta

        // FIXME: Move to library with meta programming
        template <typename F, typename T, std::size_t N, std::size_t... Idx>
        void for_all_do_impl( F f, const std::array<T, N>& a, std::index_sequence<Idx...> )
        {
            ( f( a[ Idx ] ), ... );
        } // meta

        // FIXME: Move to library with meta programming
        template <typename F, typename T, std::size_t N> void for_all_do( F f, const std::array<T, N>& a )
        {
            for_all_do_impl( f, a, std::make_index_sequence<N>{} );
        } // meta

        /** @brief Writes the content of the buffer aBuf to the output-stream using CSV.
         */
        inline void writeBufToStream( )
        {
            for_all_do( [&]( auto& rTpl ) { prtTpltoStream( ofStream, rTpl ); }, aBuf );
        } // method

        /** @brief Creates a filename for the CSV-file to be uploaded using the directory given by the backend.
         */
        fs::path uploadFileName( const std::string& rsConId, // id of the SQL connection
                                 const std::string& rsExtension ) // possible additive extension
        {
            auto pUploadDir = pHostTblDB->getDataUploadDir( );

            return pUploadDir /
                   fs::path( std::string( "upload_" ) + sHostTblName + "_con_" + rsConId + rsExtension + ".csv" );
        } // method

        // Private attributes of SQLFileBulkInserter
        std::string sHostTblName; // name of host table
        std::shared_ptr<DBCon> pHostTblDB; // database connection of the host table /--
        fs::path sCSVFileName; // name of the file for CSV output
        std::ofstream ofStream; // pointer to output stream
        bool bStreamIsVoid; // boolean state; if the stream is void, it has to be opened first
        int64_t uiReleaseThreshold; // if the file comprises uiReleaseThreshold rows, automatically release it

      public:
        SQLFileBulkInserter( SQLTable& rxHost, int64_t uiReleaseThreshold,
                             const std::string& rsExtension )
            : uiInsPos( 0 ), // initialize buffer position counter
              uiInsCnt( 0 ), // initialize counter for total number of insertions
              sHostTblName( rxHost.getTableName( ) ),
              pHostTblDB( rxHost.pDB ),
              sCSVFileName( uploadFileName( rxHost.pDB->sConId, rsExtension ) ), // filename used for CSV upload
              bStreamIsVoid( true ), // Boolean state that
              uiReleaseThreshold( uiReleaseThreshold )
        {
            std::cout << "COL NAMES:";
            for( auto rxName : rxHost.tableColNames( ) )
                std::cout << rxName << " ";
            std::cout << std::endl;
            // DEBUG: std::cout << "SQLFileBulkInserter filename:" << this->sCSVFileName.generic_string( ) <<
            // std::endl;
        } // constructor

        /** @brief Inserts a row into the table belonging to the bulk inserter.
         */
        inline void insert( const InsTypes... args )
        {
            if( bStreamIsVoid ) // 'bRelease == true' implies that the stream (file) must be closed.
            {
                // Open the stream and get any existing content purged.
                // DEBUG: std::cout << sCSVFileName << std::endl;
                ofStream.open( sCSVFileName, std::ofstream::out | std::ofstream::trunc );
                // Check whether something went wrong
                if( !ofStream )
                    throw std::runtime_error( "SQLFileBulkInserter - Opening file: " + sCSVFileName.string( ) +
                                              " failed due to reason:\n" + strerror( errno ) );
                bStreamIsVoid = false;
            } // if
            aBuf[ uiInsPos ] = std::tuple<InsTypes...>( args... ); // TO DO: Check for std::move( args )...
            if( ++uiInsPos == BUF_SIZE )
            {
                writeBufToStream( ); // write buffer to file
                uiInsPos = 0;
            } // if

            // Release if file limit reached
            if( ++uiInsCnt % ( static_cast<const size_t>( uiReleaseThreshold ) ) == 0 )
                this->release( );
        } // method

        /** @brief: Actually insert the values in the file to the table.
         *  With MySQL this results in a LOAD statement.
         */
        inline void release( bool bRelDueToDestruct = false )
        {
            if( !bStreamIsVoid )
            {
                this->flush( ); // flush remaining buffer rows
                // IMPORTANT: Status setting must come first, because fillTableByFile can throw an exception
                // (in this case, the release() call of the destructor should not go to the true branch!)
                bStreamIsVoid = true;
                // Close the stream so that it can be opened by the reader ...
                ofStream.close( );
                pHostTblDB->fillTableByFile( sCSVFileName, this->sHostTblName ); // blocks until finished
            } // if
            else if( !bRelDueToDestruct )
                throw std::runtime_error( "FileBulkInserter: Your tried to release a void stream." );
        } // method

        /** @brief Flush a non-full buffer to the file.
         */
        inline void flush( )
        {
            for( size_t uiCount = 0; uiCount < uiInsPos; uiCount++ )
                // Add the current tuple the file
                prtTpltoStream( ofStream, aBuf[ uiCount ] );

            uiInsPos = 0; // reset insert position counter
        } // method

        /** @brief Destructor flushes the buffer while dropping exceptions.
         */
        ~SQLFileBulkInserter( )
        {
            // Throwing an exception in a destructor results in undefined behavior. Therefore, we swallow these
            // exceptions and report them via std:cerr.
            doNoExcept(
                [this] {
                    this->release( true ); // true indicates that the release call is done by the destructor
                    // Delete the CSV-file if it exists.
                    if( fs::exists( this->sCSVFileName ) )
                        fs::remove( this->sCSVFileName );
                },
                "Destructor of SQLFileBulkInserter threw exception." );
        } // destructor
    }; // class (SQLFileBulkInserter)

  protected: // class SQLTable
    /** @brief get the SQL-type for each table column via the C++-type.
     *  ( C++ types are automatically translated to SQL types. )
     *  @param - Translator The delivered class that must define a C++ to SQL type translation
     *  Structured JSON dump of table def: std::cout << std::setw( 2 ) << jTableDef << std::endl;
     */
    template <class Translator, typename... ColumnTypes> struct ColTypeTranslator
    {
        // column iteration: recursive case
        template <size_t N, typename TypeCurrCol, typename... TypesRemainCols> struct TranslateTypes
        {
            static void iterate( std::vector<std::string>& rvDBColTypes )
            {
                // With respect to the Translator::template:
                // https://stackoverflow.com/questions/2105901/how-to-fix-expected-primary-expression-before-error-in-c-template-code
                rvDBColTypes.push_back( Translator::template getSQLTypeName<TypeCurrCol>( ) );
                // recursive call
                TranslateTypes<N - 1, TypesRemainCols...>::iterate( rvDBColTypes );
            } // method
        }; // struct

        // column iteration: base case (via partial specialization)
        template <typename TypeCurrCol> struct TranslateTypes<1, TypeCurrCol>
        {
            static void iterate( std::vector<std::string>& rvDBColTypes )
            {
                rvDBColTypes.push_back( Translator::template getSQLTypeName<TypeCurrCol>( ) );
            } // method
        }; // struct

        /** @brief Translate C++ types to SQL types via the Functor */
        static inline std::vector<std::string> getSQLColTypes( )
        {
            std::vector<std::string> vColsSQLTypeNames;
            TranslateTypes<sizeof...( ColumnTypes ), ColumnTypes...>::iterate( vColsSQLTypeNames );
            return vColsSQLTypeNames;
        } // method

        /** @brief Dump the SQL column types to cout. */
        static void dumpSQLColTypes( )
        {
            auto vColNames = getSQLColTypes( );
            for( auto s : vColNames )
                std::cout << s << " ";
            std::cout << std::endl;
        } // method
    }; // outer struct

    // Protected attributes of class SQLTable
    // typedef std::remove_pointer<DBCon>::type DBConBaseType;
    std::shared_ptr<DBCon> pDB; // database connection (templated) //--
    const json jTableDef; // json table definition (jTableDef must stay due to references)
    const json& rjTableCols; // reference into jTableDef
    const std::vector<std::string> vSQLColumnTypes; // SQL types of the table columns as text
    bool bDropOnDestr; // if true, table gets dropped at object destruction
    std::unique_ptr<typename DBCon::DBImplForward::PreparedStmt> pInsertStmt; // pointer to prepared insert statement

#if 0
    const std::set<size_t> xAutoNullCols; // positions of all columns receiving auto NULL insertion
#endif

    /** @brief Make the text of an appropriate SQL table creation statement.
     *  TO DO: Improved parsing and error reporting.
     */
    std::string makeTableCreateStmt( )
    {
        std::string sStmt = "CREATE TABLE " + getTableName( ) + " (";
        // deal with regular columns
        for( size_t iItr = 0; iItr < this->rjTableCols.size( ); )
        {
            auto& rjCol = this->rjTableCols[ iItr ]; // current column in jTableDef
            sStmt
                .append( rjCol[ COLUMN_NAME ] ) // column name
                .append( " " )
                .append( vSQLColumnTypes[ iItr ] ); // column type
            if( rjCol.count( CONSTRAINTS ) )
                // CONSTRAINTS must be plain text describing all constraints
                sStmt.append( " " ).append( rjCol[ CONSTRAINTS ] );
            // insert separating comma
            if( ++iItr < this->rjTableCols.size( ) )
                sStmt.append( ", " );
        } // for
        // deal with generated columns if there are any defined
        if( this->jTableDef.count( GENERATED_COLUMNS ) )
            for( size_t iItr = 0; iItr < this->jTableDef[ GENERATED_COLUMNS ].size( ); iItr++ )
            {
                // get the json column enty
                auto& rjCol = this->jTableDef[ GENERATED_COLUMNS ][ iItr ]; // current column in jTableDef

                // check that each colum is complete
                if( rjCol.count( COLUMN_NAME ) == 0 )
                    throw std::runtime_error(
                        std::string( "COLUMN_NAME is missing for a GENERATED_COLUMNS in the table " )
                            .append( getTableName( ) ) );
                if( rjCol.count( TYPE ) == 0 )
                    throw std::runtime_error( std::string( "TYPE is missing for the generated column: " )
                                                  .append( rjCol[ COLUMN_NAME ] )
                                                  .append( " in the table " )
                                                  .append( getTableName( ) ) );
                if( rjCol.count( AS ) == 0 )
                    throw std::runtime_error( std::string( "AS is missing for the generated column: " )
                                                  .append( rjCol[ COLUMN_NAME ] )
                                                  .append( " in the table " )
                                                  .append( getTableName( ) ) );
                // append to the SQL statemnt
                sStmt
                    // insert separating comma. since the last regular colum does not have a comma we can always saveley
                    // insert a seperating comma
                    .append( ", " )
                    .append( rjCol[ COLUMN_NAME ] ) // column name
                    .append( " " )
                    .append( rjCol[ TYPE ] ) // column type
                    /* Generated columns are columns that are computed from other columns
                     * InnoDb supports indices on VIRTUAL generated columns, so there is no need to make use of
                     * the STORED variant. -> actually it doesn't...
                     * The user is supposed to supply an expression that computes the generated column
                     * see:
                     * https://dev.mysql.com/doc/refman/5.7/en/create-table-generated-columns.html
                     * https://dev.mysql.com/doc/refman/5.7/en/create-table-secondary-indexes.html
                     */
                    .append( " AS ( " )
                    .append( rjCol[ AS ] )
                    .append( " ) " )
                    .append( rjCol.count( STORED ) == 0 || rjCol[ STORED ] ? "STORED" : "VIRTUAL" );
                // generated columns can have constraints as well
                if( rjCol.count( CONSTRAINTS ) )
                    // CONSTRAINTS must be plain text describing all constraints
                    sStmt.append( " " ).append( rjCol[ CONSTRAINTS ] );
            } // for

        // Primary Key (for composite form of primary key)
        if( jTableDef.count( PRIMARY_KEY ) )
            sStmt.append( ", PRIMARY KEY (" ).append( jTableDef[ PRIMARY_KEY ] ).append( ")" );

        // Foreign Key definitions
        for( auto& rjItem : jTableDef.items( ) )
            if( rjItem.key( ) == FOREIGN_KEY )
            { // Found definition
                auto& rjVal = rjItem.value( );
                if( rjVal.count( COLUMN_NAME ) != 1 || rjVal.count( REFERENCES ) != 1 )
                    throw std::runtime_error( "JSON table definition for table: " + getTableName( ) +
                                              "\nFOREIGN_KEY definition does not comprise exactly "
                                              "one COLUMN_NAME item as well as REFERENCE item." );
                sStmt.append( ", FOREIGN KEY (" )
                    .append( rjVal[ COLUMN_NAME ] )
                    .append( ")" )
                    .append( " REFERENCES " )
                    .append( rjVal[ REFERENCES ] );
            } // if
        // Close statement
        sStmt.append( ")" );
#ifdef SQL_VERBOSE
        std::cout << "SQL table creation statement: " << sStmt << std::endl;
#endif
        return sStmt;
    } // method

    /** @brief Implementation part of getValuesStmt */
    template <std::size_t... Is> void getValueStmtImpl( std::string& sStmtText, std::index_sequence<Is...> ) const
    {
        ( ( sStmtText.append( ( Is == 0 ? "" : "," ) )
                .append( DBCon::DBImplForward::TypeTranslator::template getPlaceholderForType<ColTypes>( ) ) ),
          ... );
    } // meta

    /** @brief Delivers the "VALUES"-art of an insertion statement */
    inline std::string getValuesStmt( ) const
    {
        std::string sStmtText( "(" );
        getValueStmtImpl( sStmtText, std::index_sequence_for<ColTypes...>{} );
        // For INSERT, REPLACE, and UPDATE, if a generated column is inserted into, replaced, or updated explicitly, the
        // only permitted value is DEFAULT.
        // this assumes that the create table statement always appends all generated columns at the end
        if( this->jTableDef.count( GENERATED_COLUMNS ) )
            for( size_t iItr = 0; iItr < this->jTableDef[ GENERATED_COLUMNS ].size( ); iItr++ )
                sStmtText.append( ", DEFAULT" );
        return sStmtText.append( ")" );
    } // method

  public:
    /** @brief Creates the text for a prepared SQL INSERT statement.
     *  uiNumberVals determines the number of values (rows) in the case of multiple row inserts.
     *  @details Syntax should work with most SQL dialects.
     */
    std::string makeInsertStmt( const size_t uiNumVals = 1 ) const
    {
        // Number of colums
        assert( this->rjTableCols.size( ) == sizeof...( ColTypes ) );

        const std::string sValuesPartOfStmt( getValuesStmt( ) );
        std::string sStmt = "INSERT INTO ";
        sStmt.append( getTableName( ) ).append( " VALUES " );

        for( size_t uiItrRow = 0; uiItrRow < uiNumVals; )
        {
            sStmt.append( sValuesPartOfStmt );
            // insert separating comma
            if( ++uiItrRow < uiNumVals )
                sStmt.append( ", " );
        } // outer for
#ifdef SQL_VERBOSE
        std::cout << "SQL insert statement: " << sStmt << std::endl;
#endif
        return sStmt;
    } // function
  protected:
    /** @brief Make the text of an appropriate SQL table deletion statement. */
    std::string makeTableDropStmt( )
    {
        return "DROP TABLE IF EXISTS " + getTableName( );
    } // method

    /** @brief Make a SQL statement for deleting all table rows. */
    std::string makeTableDeleteAllRowsStmt( )
    {
        return "DELETE FROM  " + getTableName( );
    } // method

    /** @brief Drop the table in DB. (Should only be done by the destructor.) */
    void drop( )
    {
        pDB->execSQL( makeTableDropStmt( ) );
    } // method

    /** @brief Get the name of the current table. */
    std::string getTableName( ) const
    {
        return jTableDef[ TABLE_NAME ];
    } // method

    /** @brief Returns the JSON definition of the table as formated string. */
    std::string dumpJsonDefToStr( )
    {
        std::stringstream ss;
        ss << "JSON table definition:\n" << std::setw( 2 ) << this->jTableDef << std::endl;
        return ss.str( );
    } // method

    /** @brief: Get the columns reference safely. */
    const json& getTableCols( )
    {
#ifdef SQL_VERBOSE
        std::cout << "JSON table definition:\n" << std::setw( 2 ) << this->jTableDef << std::endl;
#endif
        if( this->jTableDef.count( TABLE_COLUMNS ) != 1 )
            throw std::runtime_error( dumpJsonDefToStr( ) +
                                      "\nJSON table definitions require exactly one TABLE_COLUMNS section." );
        return this->jTableDef[ TABLE_COLUMNS ];
    } // method

    /** @brief: Returns a vector comprising the names of table columns. */
    std::vector<std::string> tableColNames( )
    {
        std::vector<std::string> vColNames;
        for( auto& rxCol : rjTableCols )
            if( rxCol.count( COLUMN_NAME ) )
                vColNames.push_back( rxCol[ COLUMN_NAME ] );
            else
                throw std::runtime_error( "JSON for table definition: Missing item COLUMN_NAME.\n" +
                                          dumpJsonDefToStr( ) );

        return vColNames;
    } // method

  public:
    typedef ColTypeTranslator<typename DBCon::TypeTranslator, ColTypes...> CollTypeTranslation;

    /** @brief Initializes the table for use.
     *  @detail If the table does not exist in the DB, the table is created first.
     */
    void init( )
    {
        // Create table if not existing in DB guarantee
        pDB->doPoolSafe( [&] {
            if( !pDB->tableExistsInDB( getTableName( ) ) )
                pDB->execSQL( makeTableCreateStmt( ) );
            else
            {
                // Verify the table for correctness
            } // else
        } );

        // Create prepared insert statement, now where the existence of the table is guaranteed
        // Improvement: Could be done on demand ...
        pInsertStmt = std::make_unique<typename DBCon::PreparedStmt>( this->pDB, makeInsertStmt( ) );

        // Parse CPP extra section
        if( jTableDef.count( CPP_EXTRA ) )
        {
            for( auto& rxCPP_Extra : jTableDef[ CPP_EXTRA ] )
            {
                if( rxCPP_Extra == "DROP ON DESTRUCTION" )
                {
                    this->bDropOnDestr = true;
                } // if
            } // for
        } // if

        // INACTIVE: // Parse SQL extra section
        // INACTIVE: if( jTableDef.count( SQL_EXTRA ) )
        // INACTIVE: {
        // INACTIVE:     for( auto& rxSQL_Extra : jTableDef[ SQL_EXTRA ] )
        // INACTIVE:     {
        // INACTIVE:     } // for
        // INACTIVE: } // if
    } // method

    /** @brief Dump the table content. */
    void dump( )
    {
        SQLQuery<DBCon, ColTypes...> xQuery( this->pDB, "SELECT * FROM " + getTableName( ) );
        xQuery.execAndForAll( []( const auto&... val ) {
            ( ( std::cout << val << ',' ), ... );
            std::cout << std::endl;
        } );
        // xQuery.applyToAllRows( []( const auto&... val ) { ( ( std::cout << val << ',' ), ... ); } );
    } // method


    /** @brief Insert the argument pack as fresh row into the table without guaranteeing type correctness.
     *  @detail This form of the insert does not guarantees that there are no problems with column types at runtime.
     *  It allows the injection of NULL values into an insertion by passing (void *)NULL at th column position,
     * where the NULL value shall occurs. Don't use insertNonSafe, if you do not need explicit NULL values in the
     * database.
     */
    template <typename... ArgTypes> inline void insertNonSafe( const ArgTypes&... args )
    {
        // static_assert( sizeof...( ArgTypes ) == sizeof...( ColTypes ) );
        auto iAffectedRows = pInsertStmt->bindAndExec( args... );
        if( iAffectedRows != 1 )
            throw std::runtime_error( "SQL DB insert. iAffectedRows != 1. (" + std::to_string( iAffectedRows ) + ")" );
    } // method

    /** @brief: Type-safely insert the argument pack as fresh row into the table. */
    inline SQLTable& insert( const ColTypes&... args )
    {
        this->insertNonSafe( args... );
        return *this;
    } // method

    // Alias for bulk inserter type
    template <size_t BUF_SIZE> using SQLBulkInserterType = SQLBulkInserter<BUF_SIZE, ColTypes...>;

    /** @brief Delivers a shared pointer to a bulk inserter for inserting table rows type-safely.
     *   TIP: Definition of a bulk inserter (don't forget the keyword template with GCC!):
     *        auto xBulkInserter = table.template getBulkInserter< size >();
     */
    template <size_t BUF_SIZE> std::shared_ptr<SQLBulkInserterType<BUF_SIZE>> getBulkInserter( )
    {
        return std::make_shared<SQLBulkInserterType<BUF_SIZE>>( *this );
    } // method

    /** @brief Alias for file bulk inserter type */
    template <size_t BUF_SIZE>
    using SQLFileBulkInserterType =
        typename SQLTable<DBCon, ColTypes...>::template SQLFileBulkInserter<BUF_SIZE, ColTypes...>;

    /** @brief Bulk inserter that uses a file for data steaming. */
    template <size_t BUF_SIZE>
    std::shared_ptr<SQLFileBulkInserterType<BUF_SIZE>> getFileBulkInserter( const std::string& pCSVFileName )
    {
        return std::make_shared<SQLFileBulkInserterType<BUF_SIZE>>( *this, pCSVFileName );
    } // method

    /** @brief: Delete all rows in the table. */
    SQLTable& deleteAllRows( )
    {
        pDB->execSQL( makeTableDeleteAllRowsStmt( ) );
        return *this;
    } // method

    /** @brief: Creates requested index if not already existing. */
    void addIndex( const json& rjIndexDef )
    {
        SQLIndexView test( this, rjIndexDef );
    } // method

    SQLTable( const SQLTable& ) = delete; // no table copies

    /* Constructor
     * Design improvement: The positions of the NULL insertions could be coded via the JSON table definition
     */
    SQLTable( std::shared_ptr<DBCon> pDB, // DB connector //--
              const json& rjTableDef // table definition in JSON
#if 0
              , const std::set<size_t>& xAutoNullCols = { } // positions of columns with NULL insertions
#endif
              )
        : pDB( pDB ),
          jTableDef( rjTableDef ), // deep copy
          rjTableCols( getTableCols( ) ),
          vSQLColumnTypes( CollTypeTranslation::getSQLColTypes( ) ),
          bDropOnDestr( false )
#if 0
          , xAutoNullCols( xAutoNullCols )
#endif
    {
        if( !( ( rjTableCols.size( ) == vSQLColumnTypes.size( ) ) &&
               ( rjTableCols.size( ) == sizeof...( ColTypes ) ) ) )
            throw std::runtime_error( "Incorrect SQL table definition. Number of C++ types does not match number of "
                                      "columns in JSON definition." );
        init( );
    }; // constructor

    /* Destructor */
    ~SQLTable( )
    {
        do_exception_safe( [&]( ) {
            if( pDB && bDropOnDestr )
                this->drop( );
        } );
    } // destructor
}; // class (SQLTable)


/** @brief Variation of the normal SQL table having an automatic primary key.
 *  Primary Key is of type int64_t and on an additional first column called ID.
 */
template <typename DBCon, typename... ColTypes>
#ifdef PRIMARY_KEY_ON_LAST_ROW
class SQLTableWithAutoPriKey : public SQLTable<DBCon, ColTypes..., int64_t>
#else
class SQLTableWithAutoPriKey : public SQLTable<DBCon, int64_t, ColTypes...>
#endif
{
    /* Inject primary key column definition into table definition */
    json inject( const json& rjTableDef )
    {
        // Create deep copy of the table definition
        json jTableDef( rjTableDef );
#ifdef PRIMARY_KEY_ON_LAST_ROW
        if( jTableDef.count( TABLE_COLUMNS ) )
            jTableDef[ TABLE_COLUMNS ].push_back( json::object( {{COLUMN_NAME, "ID"},
                                                                 { CONSTRAINTS,
                                                                   "NOT NULL AUTO_INCREMENT PRIMARY KEY" }} ) );
#else
        // Primary key on first row:
        if( jTableDef.count( TABLE_COLUMNS ) )
            jTableDef[ TABLE_COLUMNS ].insert(
                jTableDef[ TABLE_COLUMNS ].begin( ),
                json::object( {{COLUMN_NAME, "id"}, {CONSTRAINTS, "NOT NULL AUTO_INCREMENT UNIQUE PRIMARY KEY"}} ) );
#endif
        return jTableDef;
    } // method

    /** @brief Insert the argument pack as fresh row into the table in a thread-safe way.
     *	Delivers the primary key of the inserted row as return value.
     *  @detail We need a lock_guard, because in the case of concurrent inserts (on DB level) we have to guarantee
     *that all AUTO_INCREMENT inserts occur serialized for getting the correct AUTO_INCREMENT value.
     */
    template <typename... ArgTypes> inline int64_t insertThreadSafe( const ArgTypes&... args )
    {
        std::lock_guard<std::mutex> xGuard( this->pDB->pGlobalInsertLock );
        // Thread-safe region (till end of method).
        SQLTable<DBCon, int64_t, ColTypes...>::insertNonSafe( nullptr, args... );
        return this->pDB->getLastAutoIncrVal( );
    } // method

  public:
    /**
     * @brief holds the columns types.
     * @details
     * since a parameter pack cannot be held in a using, we hold the type pack<ColTypes...>.
     * pack is a completely empty struct.
     * @see get_inserter_container_module.h, there this is used to extract the types of all column from a tabletype that
     * is given as a template parameter.
     * columnTypes -> ForwardedColTypes
     */
    using columnTypes = pack<ColTypes...>;

    SQLTableWithAutoPriKey( const SQLTableWithAutoPriKey& ) = delete; // no table copies

    /* Constructor */
    SQLTableWithAutoPriKey( std::shared_ptr<DBCon> pDB, const json& rjTableDef ) //--
#ifdef PRIMARY_KEY_ON_LAST_ROW
        : SQLTable<DBCon, ColTypes..., int64_t>( pDB, inject( rjTableDef ),
                                                 std::set<size_t>( {rjTableDef[ TABLE_COLUMNS ].size( )} )
#else
        : SQLTable<DBCon, int64_t, ColTypes...>( pDB, inject( rjTableDef ) /* , std::set<size_t>( { 0 } ) */
#endif
          )
    {} // constructor

    /** @brief Inserts the argument pack as fresh row into the table without guaranteeing type correctness.
     *	Delivers the primary key of the inserted row as return value.
     *  @detail It allows the injection of NULL values into an insertion by passing nullptr at the column position,
     *   where the NULL value shall occurs.
     *   Don't use insertNonSafe, if you do not need explicit NULL values in the database.
     */
    template <typename... ArgTypes> inline int64_t insertNonSafe( const ArgTypes&... args )
    {
        return this->insertThreadSafe( args... );
    } // method

    /** @brief Inserts a fresh row into the table type-safely. (type-safe variant of insertNonSafe)
     * Delivers the primary key of the inserted row as return value.
     * FIXME: std::forward<ArgTypes>( args )...
     */
    inline int64_t insert( const ColTypes&... args )
    {
        return this->insertThreadSafe( args... );
    } // method

    /** @brief Alias for bulk inserter type for tables with primary key */
    template <size_t BUF_SIZE>
    using SQLBulkInserterType =
        typename SQLTable<DBCon, int64_t, ColTypes...>::template SQLBulkInserter<BUF_SIZE, std::nullptr_t, ColTypes...>;

    /** @brief Delivers a shared pointer to a bulk inserter for inserting table rows type-safely.
     *   TIP: Definition of a bulk inserter (don't forget the keyword template with GCC!):
     *        auto xBulkInserter = table.template getBulkInserter< size >();
     */
    template <size_t BUF_SIZE> std::shared_ptr<SQLBulkInserterType<BUF_SIZE>> getBulkInserter( )
    {
        return std::make_shared<SQLBulkInserterType<BUF_SIZE>>( *this );
    } // method

    /** @brief Alias for file bulk inserter type for tables with primary key */
    template <size_t BUF_SIZE>
    using SQLFileBulkInserterType =
        typename SQLTable<DBCon, int64_t, ColTypes...>::template SQLFileBulkInserter<BUF_SIZE, std::nullptr_t,
                                                                                     ColTypes...>;

    /** @brief Bulk inserter that uses a file for data steaming. */
    template <size_t BUF_SIZE>
    std::shared_ptr<SQLFileBulkInserterType<BUF_SIZE>> getFileBulkInserter( int64_t uiReleaseThreshold = 1000000,
                                                                            const std::string& srExtension = "" )
    {
        return std::make_shared<SQLFileBulkInserterType<BUF_SIZE>>( *this, uiReleaseThreshold, srExtension );
    } // method
}; // class
