/* Authors: Arne Kutzner and Markus Schmidt
 * Created: Jan. 2020
 * MIT License
 * @file common.h
 * @brief General API for SQL-like databases. Use this API in combination with a connector to work with SQL databases.
 */

#pragma once

// #define SQL_VERBOSE // define me if required
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

#include <json.hpp>
// json support should have been included already by the connector ...
using nlohmann::json;

/** @brief Returns a string that contains info regarding types and values of an argument pack.
 *  Primarily for debugging purposes...
 */
template <typename... ArgTypes> std::string dumpArgPack( ArgTypes&&... args )
{
    std::stringstream xInfStream;
    ( ( xInfStream << typeid( args ).name( ) << ":" << args << "  " ), ... );
    return xInfStream.str( );
}

/** @brief An instance represents a connection to DB-system like MySQL etc.
 *   DBImpl must be a class that implements a database API.
 *  @details A single connection is not threadsafe.
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
    }; // class

    std::shared_ptr<bool> pTombStone; // So long the database connection is alive the stone says false.
                                      // As soon as the DB is gone, the stone says true.

  public:
    typedef std::unique_ptr<GuardedTransaction> uniqueGuardedTrxnType;
    typedef std::shared_ptr<GuardedTransaction> sharedGuardedTrxnType;

    SQLDB( const SQLDB& ) = delete; // no DB Connection copies

    // For serializing all table insertions with respect to the current DB-connection.
    // Necessary for getting the correct last insertion ID in case of an auto-increment column. (e.g. with MySQL)
    std::mutex pGlobalInsertLock;

    SQLDB( ) : DBImpl( ), pTombStone( std::make_shared<bool>( false ) )
    {}

    SQLDB( const std::string& rsDBName )
    {
        throw std::runtime_error( "DB-Construction with name not implemented jet." );
    } // constructor

    ~SQLDB( )
    {
        // Inform about the "death" of the database connection.
        *pTombStone = true;
        // Now the destructor of the DB implementor is called ...
    }; // destructor

    int64_t lastRowId( )
    {
        throw std::runtime_error( "lastRowId of SQLDB is not implemented yet." );
    } // method

    /** @brief Directly execute the statement in DB.
     */
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

    std::unique_ptr<GuardedTransaction> sharedGuardedTrxn( )
    {
        return std::make_unique<GuardedTransaction>( *this, pTombStone );
    } // method

    // FIXME: Not safe now in the case of several connections.
    bool indexExists( const std::string& rsTblName, const std::string& rsIdxName )
    {
        return DBImpl::indexExistsInDB( rsTblName, rsIdxName );
    } // method
}; // SQLDB


template <typename DBCon> class SQLStatement
{
  protected:
    std::shared_ptr<SQLDB<DBCon>> pDB; // Database connection
    std::unique_ptr<typename DBCon::PreparedStmt> pStmt; // pointer to insert statement

  public:
    SQLStatement( const SQLStatement& ) = delete; // no statement copies

    SQLStatement( std::shared_ptr<SQLDB<DBCon>> pDB, const std::string& rsStmtText )
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
    std::shared_ptr<SQLDB<DBCon>> pDB; // Database connection
    std::unique_ptr<typename DBCon::template PreparedQuery<ColTypes...>> pQuery; // pointer to insert statement
    const std::string sStmtText; // backup of the query statement text. (used for verbose error reporting)
    /* Deprecated. Execute and bind the query args., but do not fetch the first row. */

    template <typename... ArgTypes> void execAndBind( ArgTypes&&... args )
    {
        // Execute statement
        pQuery->bindAndExec( std::forward<ArgTypes>( args )... );
        // Bind for later fetching
        pQuery->bindResult( );
    } // method

    template <typename F> void forAllRowsDo( F&& func )
    {
        // Fetch all rows and apply func to all of them
        while( pQuery->fetchNextRow( ) )
            // unpack the tuple and call func (C++17 fold expression)
            STD_APPLY( func, pQuery->tCellValues );
    } // method

  public:
    SQLQuery( std::shared_ptr<SQLDB<DBCon>> pDB, const std::string& rsStmtText )
        : pDB( pDB ),
          pQuery( std::make_unique<typename DBCon::template PreparedQuery<ColTypes...>>( this->pDB, rsStmtText ) ),
          sStmtText( rsStmtText )
    {} // constructor

    template <typename... ArgTypes> bool execAndFetch( ArgTypes&&... args )
    {
        // Execute the query
        this->execAndBind( std::forward<ArgTypes>( args )... );
        // Fetch first row for getting correct EOF value.
        return this->next( );
    } // method

    /** @brief Execute the query passing args and call the function func for all rows of the query outcome
     */
    template <typename F, typename... ArgTypes> void execAndForAll( F&& func, ArgTypes&&... args )
    {
        this->execAndBind( std::forward<ArgTypes>( args )... );
        this->forAllRowsDo( std::forward<F>( func ) );
    } // method

    /** @brief Execute the query and get the n-th cell of the first row of a query.
     */
    template <int COL_NUM, typename... ArgTypes>
    typename std::tuple_element<COL_NUM, std::tuple<ColTypes...>>::type execAndGetNthCell( ArgTypes&&... args )
    {
        // Execute query, bind, and fetch first row
        this->execAndBind( std::forward<ArgTypes>( args )... );
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

    /* Delivers a reference to a tuple holding the data of the current row.
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
     *         This kind of eof works in a "look back" fashion; it tells if the previous fetch failed or not.
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
        this->execAndBind( std::forward<ArgTypes>( args )... );
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
        this->execAndBind( std::forward<ArgTypes>( args )... );
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
const std::string CPP_EXTRA = "CPP_EXTRA";
const std::string SQL_EXTRA = "SQL_EXTRA";
const std::string PRIMARY_KEY = "PRIMARY_KEY";
const std::string FOREIGN_KEY = "FOREIGN_KEY";
const std::string REFERENCES = "REFERENCES";

// Constants for index definitions via json
const std::string INDEX_NAME = "INDEX_NAME";
const std::string INDEX_COLUMNS = "INDEX_COLUMNS";
const std::string WHERE = "WHERE";

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
         * Tuple elements are moved for efficiency, so the array and its elements should not be moved after catenation.
         */
        template <typename T, std::size_t N, typename Indices = std::make_index_sequence<N>>
        auto cat_arr_of_tpl( const std::array<T, N>& a )
        {
            return cat_arr_of_tpl_impl( a, Indices{ } );
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
            // This statement creates trouble for some GCC compilers (template instantiation depth exceeds maximum) ...
            STD_APPLY( [ & ]( auto&... args ) { pBulkInsertStmt->bindAndExec( args... ); }, tCatTpl );
        } // method
#else // Iterating approach for parameter binding with bulk inserts. (Lacks beauty, but better manageable for compilers)

        /** @brief Unpack the tuple to an argument pack and bind using this pack.
         *  (Part of bulkInsert)
         */
        template <int OFFSET, typename... Args> void doSingleBind( const std::tuple<Args...>& rTpl )
        {
            STD_APPLY( [ & ]( auto&... args ) //
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
            forAllDoBindImpl( a, std::make_index_sequence<N>{ } );
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
        SQLBulkInserter( SQLTable& rxHost )
            : uiInsPos( 0 ),
              pBulkInsertStmt( // compile bulk insert stmt
                  std::make_unique<typename DBCon::PreparedStmt>( rxHost.pDB, rxHost.makeInsertStmt( BUF_SIZE ) ) ),
              pSingleInsertStmt( // compile single insert stmt
                  std::make_unique<typename DBCon::PreparedStmt>( rxHost.pDB, rxHost.makeInsertStmt( 1 ) ) )
        {
            static_assert( sizeof...( InsTypes ) == sizeof...( ColTypes ) );
        } // constructor

        /** @brief Inserts a row into the table via a bulk-insert approach.
         */
        inline void insert( const InsTypes... args )
        {
            aBuf[ uiInsPos ] = std::tuple<InsTypes...>( args... ); // TO DO: Check for std::move( args )...
            if( ++uiInsPos == BUF_SIZE )
            {
                // bulkInsert( ) does not inspects uiInsPos anymore. However, it can throw an exceptions.
                // In order to avoid the crash scenario described in flush(),  uiInsPos must be set zero before calling
                // bulkInsert().
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
                STD_APPLY( [ & ]( auto&... args ) { pSingleInsertStmt->bindAndExec( args... ); }, aBuf[ uiCount ] );
        } // method

        /** @brief Destructor flushes the buffer.
         */
        ~SQLBulkInserter( )
        {
            this->flush( ); // Due to the design, we do not need to swallow exceptions here (see flush).
        } // destructor
    }; // class

    template <size_t BUF_SIZE,
              typename... InsTypes> // InsTypes - either equal to corresponding type in ColTypes or std::nullptr_t
    class SQLFileBulkInserter
    {
      private:
        std::array<std::tuple<InsTypes...>, BUF_SIZE> aBuf; // buffer of values inserted by a single op.
        size_t uiInsPos; // index of the next insert into the array
        size_t uiInsCnt; // counter for all inserts so far

        /** @brief Serializes a value for CSV representation.
         *  NULL values are represented by '\N'.
         */
        template <typename Type> std::string csvPrint( Type& rVal )
        {
            return std::to_string( rVal );
        }; // method

        std::string csvPrint( const std::string& rVal )
        {
            return "'" + rVal + "'";
        }; // method

        std::string csvPrint( const std::nullptr_t& rVal )
        {
            return "\\N";
        }; // method

        // Tuple printer for ostream.
        // Found at: https://en.cppreference.com/w/cpp/utility/integer_sequence
        template <class Ch, class Tr, class Tuple, std::size_t... Is>
        void prtTpltoStreamImpl( std::basic_ostream<Ch, Tr>& os, const Tuple& t, std::index_sequence<Is...> )
        {
            // Requires C++17 folded expression support ...
            ( ( os << ( Is == 0 ? "" : "\t" ) << csvPrint( std::get<Is>( t ) ) ), ... );
        } // meta

        template <class Ch, class Tr, class... Args>
        void prtTpltoStream( std::basic_ostream<Ch, Tr>& os, const std::tuple<Args...>& t )
        {
            prtTpltoStreamImpl( os, t, std::index_sequence_for<Args...>{ } );
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
            for_all_do_impl( f, a, std::make_index_sequence<N>{ } );
        } // meta

        /** @brief Writes the content of the buffer aBuf to the output-stream using CSV.
         */
        inline void writeBufToStream( )
        {
            for_all_do( [ & ]( auto& rTpl ) { prtTpltoStream( ofStream, rTpl ); }, aBuf );
        } // method

        fs::path uploadFileName( const std::string& rsExtension )
        {
            std::string sSecureFilePrivPath = pHostTblDB->getVarSecureFilePriv( );
			if (sSecureFilePrivPath.empty())
				// sSecureFilePrivPath = "/MAdata/tmp";
                sSecureFilePrivPath = "/tmp";

            return fs::path(sSecureFilePrivPath) /
                   fs::path( std::string( "upload_" ) + sHostTblName + rsExtension + ".csv" );
        } // method

        std::string sHostTblName; // name of host table
        std::shared_ptr<SQLDB<DBCon>> pHostTblDB; // database connection of the host table
        fs::path sCSVFileName; // name of the file for CSV output
        std::ofstream ofStream; // pointer to output stream
        bool bStreamIsVoid; // boolean state. If stream is void, it has to be opened first

      public:
        SQLFileBulkInserter( SQLTable& rxHost, const std::string& sCSVFileName )
            : sHostTblName( rxHost.getTableName( ) ),
              pHostTblDB( rxHost.pDB ),
              sCSVFileName( uploadFileName( sCSVFileName ) ), // must become a JSON-parameter
              bStreamIsVoid( true ) // Boolean state that
        {
            // Quick hack:
            std::cout << "FILE NAME:" << this->sCSVFileName.generic_string( ) << std::endl;
        } // constructor

        /** @brief Inserts a row into the table belonging to the bulk inserter.
         */
        inline void insert( const InsTypes... args )
        {
            if( bStreamIsVoid ) // 'bRelease == true' implies that the stream (file) must be closed.
            {
                // Open the stream and get any existing content purged.
                ofStream.open( sCSVFileName, std::ofstream::out | std::ofstream::trunc );
                bStreamIsVoid = false;
            } // if
            aBuf[ uiInsPos ] = std::tuple<InsTypes...>( args... ); // TO DO: Check for std::move( args )...
            if( ++uiInsPos == BUF_SIZE )
            {
                writeBufToStream( ); // write buffer to file
                uiInsPos = 0;
            } // if

            // Release if file limit reached
            if( ++uiInsCnt % 100000 == 0 )
                this->release( );
        } // method

        /** @brief: Actually insert the values in the file to the table.
         *  With MySQL this results in a LOAD statement.
         */
        inline void release( bool bRelDoToDestruct = false )
        {
            if( !bStreamIsVoid )
            {
                this->flush( ); // flush before release
                // IMPORTANT: Status setting must come first, because fillTableByFile can throw an exception
                // (in this case, the release() call of the destructor must not go to the true branch!)
                bStreamIsVoid = true;
                // Close the stream so that it can be opened by the reader
                ofStream.close( );

                try
                {
                    pHostTblDB->fillTableByFile( sCSVFileName, this->sHostTblName ); // blocks until finished
                }
                catch( std::exception& e )
                {
                    std::cout << "EX:" << e.what( ) << std::endl;
                }
                // Close the stream so that it can be reopened ...

            } // if
            else if( !bRelDoToDestruct )
                throw std::runtime_error( "FileBulkInserter: Your tried to release a void stream." );
        } // method

        /** @brief Flush a non-full buffer to the database.
         */
        inline void flush( )
        {
            for( size_t uiCount = 0; uiCount < uiInsPos; uiCount++ )
                // Add the current tuple the file
                prtTpltoStream( ofStream, aBuf[ uiCount ] );

            uiInsPos = 0; // reset insert position counter
        } // method

        /** @brief Destructor flushes the buffer
         */
        ~SQLFileBulkInserter( )
        {
            this->release( true ); // true indicates that the release call is done by the destructor
            // Delete the CSV-file if it exists.
            // REINSERT: if( fs::exists( this->sCSVFileName ) )
            // REINSERT:     fs::remove( this->sCSVFileName );
        } // destructor
    }; // class (SQLFileBulkInserter)

  protected:
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

    // Protected attributes
    std::shared_ptr<SQLDB<DBCon>> pDB; // database connection (templated)
    const json jTableDef; // json table definition (jTableDef must stay due to references)
    const json& rjTableCols; // reference into jTableDef
    const std::vector<std::string> vSQLColumnTypes; // SQL types of the table columns as text
    bool bDropOnDestr; // if true, table gets dropped at object destruction
    std::unique_ptr<typename DBCon::PreparedStmt> pInsertStmt; // pointer to prepared insert statement

#if 0
    const std::set<size_t> xAutoNullCols; // positions of all columns receiving auto NULL insertion
#endif

    /** @brief Make the text of an appropriate SQL table creation statement.
     *  TO DO: Improved parsing and error reporting.
     */
    std::string makeTableCreateStmt( )
    {
        std::string sStmt = "CREATE TABLE " + getTableName( ) + " (";
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
  public:
    /** @brief Creates the text for a prepared SQL INSERT statement.
     *  uiNumberVals determines the number of values (rows) in the case of multiple row inserts.
     *  @details Syntax should work with most SQL dialects.
     */
    std::string makeInsertStmt( const size_t uiNumVals = 1 )
    {
        std::string sStmt = "INSERT INTO ";
        sStmt.append( getTableName( ) ).append( " VALUES " );

        for( size_t uiItrRow = 0; uiItrRow < uiNumVals; )
        {
            sStmt.append( "(" );
            // Comma separated list
            for( size_t uiItr = 0; uiItr < this->rjTableCols.size( ); )
            {
#if 0
            sStmt.append( this->xAutoNullCols.find( uiItr ) != this->xAutoNullCols.end( ) ? "NULL" : "?" );
#else
                sStmt.append( "?" );
#endif
                // insert separating comma
                if( ++uiItr < this->rjTableCols.size( ) )
                    sStmt.append( ", " );
            } // inner for
            sStmt.append( ")" );
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

    std::string dumpJsonDefToStr( )
    {
        std::stringstream ss;
        ss << "JSON table definition:\n" << std::setw( 2 ) << this->jTableDef << std::endl;
        return ss.str( );
    } // method

    /** @brief: Get the columns reference safely.
     */
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

  public:
    typedef ColTypeTranslator<typename DBCon::TypeTranslator, ColTypes...> CollTypeTranslation;

    /** @brief Initializes the table for use.
     *  @details If the table does not exist in the DB, the table is created first.
     */
    void init( )
    {
        // Create table if not existing in DB
        if( !pDB->tableExistsInDB( getTableName( ) ) )
            pDB->execSQL( makeTableCreateStmt( ) );
        else
        {
            // Verify the table for correctness
        }

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
     *  @details This form of the insert does not guarantees that there are no problems with column types at runtime.
     *  It allows the injection of NULL values into an insertion by passing (void *)NULL at th column position, where
     *  the NULL value shall occurs.
     *  Don't use insertNonSafe, if you do not need explicit NULL values in the database.
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
    SQLTable( std::shared_ptr<SQLDB<DBCon>> pDB, // DB connector
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
        if( rjTableCols.size( ) != vSQLColumnTypes.size( ) )
            throw std::runtime_error(
                "Incorrect SQL table definition. Number of C++ types does not match number of columns." );
        init( );
    }; // constructor

    /* Destructor */
    ~SQLTable( )
    {
        do_exception_safe( [ & ]( ) {
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
            jTableDef[ TABLE_COLUMNS ].push_back( json::object( { { COLUMN_NAME, "ID" },
                                                                  { CONSTRAINTS,
                                                                    "NOT NULL AUTO_INCREMENT PRIMARY KEY" } } ) );
#else
        // Primary key on first row:
        if( jTableDef.count( TABLE_COLUMNS ) )
            jTableDef[ TABLE_COLUMNS ].insert(
                jTableDef[ TABLE_COLUMNS ].begin( ),
                json::object(
                    { { COLUMN_NAME, "id" }, { CONSTRAINTS, "NOT NULL AUTO_INCREMENT UNIQUE PRIMARY KEY" } } ) );
#endif
        return jTableDef;
    } // method

    /** @brief Insert the argument pack as fresh row into the table in a thread-safe way.
     *	Delivers the primary key of the inserted row as return value.
     *  @details We need a lock_guard, because in the case of concurrent inserts (on DB level) we have to guarantee that
     *   all AUTO_INCREMENT inserts occur serialized for getting the correct AUTO_INCREMENT value.
     */
    template <typename... ArgTypes> inline int64_t insertThreadSafe( const ArgTypes&... args )
    {
        std::lock_guard<std::mutex> xGuard( this->pDB->pGlobalInsertLock );
        // Thread-safe region (till end of method).
        SQLTable<DBCon, int64_t, ColTypes...>::insertNonSafe( nullptr, args... );
        return this->pDB->getLastAutoIncrVal( );
    } // method

  public:
    SQLTableWithAutoPriKey( const SQLTableWithAutoPriKey& ) = delete; // no table copies

    /* Constructor */
    SQLTableWithAutoPriKey( std::shared_ptr<SQLDB<DBCon>> pDB, const json& rjTableDef )
#ifdef PRIMARY_KEY_ON_LAST_ROW
        : SQLTable<DBCon, ColTypes..., int64_t>( pDB, inject( rjTableDef ),
                                                 std::set<size_t>( { rjTableDef[ TABLE_COLUMNS ].size( ) } )
#else
        : SQLTable<DBCon, int64_t, ColTypes...>( pDB, inject( rjTableDef ) /* , std::set<size_t>( { 0 } ) */
#endif
          )
    {} // constructor

    /** @brief Inserts the argument pack as fresh row into the table without guaranteeing type correctness.
     *	Delivers the primary key of the inserted row as return value.
     *  @details It allows the injection of NULL values into an insertion by passing nullptr at the column position,
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

    // Alias for bulk inserter type
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

    // Alias for bulk inserter type
    template <size_t BUF_SIZE>
    using SQLFileBulkInserterType =
        typename SQLTable<DBCon, int64_t, ColTypes...>::template SQLFileBulkInserter<BUF_SIZE, std::nullptr_t,
                                                                                     ColTypes...>;

    template <size_t BUF_SIZE>
    std::shared_ptr<SQLFileBulkInserterType<BUF_SIZE>> getFileBulkInserter( const std::string& pCSVFileName )
    {
        return std::make_shared<SQLFileBulkInserterType<BUF_SIZE>>( *this, pCSVFileName );
    } // method
}; // class
