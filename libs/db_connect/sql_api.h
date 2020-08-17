/* Authors: Arne Kutzner and Markus Schmidt
 * Created: Jan. 2020
 * MIT License
 * @file sql_api.h
 * @brief General API for SQL-like databases. Use this API in combination with a connector to work with SQL databases.
 */

#pragma once
// #define SQL_VERBOSE // define me if required
// #define POSTGRESQL

#ifdef _MSC_VER
#pragma warning( disable : 4996 ) // suppress warnings regarding the use of strerror
#endif

#include "fort.hpp"
#include "util/support.h"
#include <atomic> // atomic fetch_add
#include <cerrno> // error management for file I/O
#include <cstring> // error management for file I/O
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <vector>

#ifndef DB_BASE_INCLUDED
#include <db_base.h>
#endif

// For getting the code working with gcc 6.x.x compiler
#if defined( __GNUC__ ) && ( __GNUC__ < 7 ) && !defined( __clang__ )
#include <experimental/tuple>
#define STD_APPLY std::experimental::apply
#else
#include <tuple>
#define STD_APPLY std::apply
#endif

#if( defined( __GNUC__ ) && ( __GNUC__ < 8 ) )
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

// Import JSON support
#include <json.hpp>
using nlohmann::json;

/**
 * @brief completely empty struct
 * @details
 * can be used to hold and then extract a parameter TypePack from a variable.
 */
template <typename... Args> struct TypePack
{};

/** @brief Required for workaround with respect to explicit specialization of a template function in a template class.
 * See: https://stackoverflow.com/questions/3052579/explicit-specialization-in-non-namespace-scope
 */
template <typename Type> struct Identity
{
    typedef Type TemplatedType;
}; // struct

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

using PriKeyDefaultType = uint64_t; // Default type of a primary key columns that are managed by the library

/** @brief Class definition for single master sync object managing concurrent access and primary key counters. */
class SQLDBGlobalSync
{
    using CntsMapType = std::map<std::string, std::shared_ptr<std::atomic<PriKeyDefaultType>>>;

    std::shared_ptr<CntsMapType> pPriKeyCntsMap; // master map that keeps associations
                                                 // between table names and primary key counters
    std::unique_ptr<std::mutex> pMasterMutex; // top level mutex of SQL API.

    // Map the counts for schema the number of connections to them
    std::unique_ptr<std::map<std::string, uint32_t>> pSchemaCntsMap;

  public:
    /** @brief Constructs and initializes the global sync object */
    SQLDBGlobalSync( )
        : pPriKeyCntsMap( std::make_shared<CntsMapType>( ) ),
          pMasterMutex( std::make_unique<std::mutex>( ) ),
          pSchemaCntsMap( std::make_unique<std::map<std::string, uint32_t>>( ) )
    {
        // DEBUG: std::cout << "Init SQLDBGlobal ..." << std::endl;
    } // constructor

    /** @brief Destructs the global sync object */
    ~SQLDBGlobalSync( )
    {
        // DEBUG: std::cout << "Finish SQLDBGlobal ..." << std::endl;
    } // destructor

    /** @brief Registers a connection for the given schema with the visor. */
    void registerSchema( const std::string& rsSchemaName ) const
    {
        // Make sure that the below code is exclusively executed by one thread.
        std::lock_guard<std::mutex> xLock( *pMasterMutex );

        auto pItr = pSchemaCntsMap->find( rsSchemaName );
        if( pItr == pSchemaCntsMap->end( ) )
            // There is no primary key entry for rsSchemaName so far. Create one ...
            pSchemaCntsMap->emplace( rsSchemaName, 1 );
        else
            // Increment counter
            ( pItr->second )++;
    } // method

    /** @brief Unregisters a connection for the given schema with the visor.
     *  Returns the number of remaining connections after unregistering.
     */
    uint32_t unregisterSchema( const std::string& rsSchemaName ) const
    {
        // Make sure that the below code is exclusively executed by one thread.
        std::lock_guard<std::mutex> xLock( *pMasterMutex );

        // Decrement the counter and return the result.
        assert( pSchemaCntsMap->at( rsSchemaName ) > 0 );
        return --( pSchemaCntsMap->at( rsSchemaName ) );
    } // method

    /** @brief Returns a shared pointer to the primary key counter for the table with name rsTblName. If there is no
     *  such pointer, it is created and initialized using uiPriKeyCnt.
     *  @details Design: All tables with an automatic primary key call this method in their constructor delivering the
     *  maximum of the id column as suggestion.
     */
    std::shared_ptr<std::atomic<PriKeyDefaultType>> getPtrToPriKeyCnt( const std::string rsTblName,
                                                                       PriKeyDefaultType uiPriKeyCnt ) const
    {
        // Make sure that the below code is exclusively executed by one thread.
        std::lock_guard<std::mutex> xLock( *pMasterMutex );

        auto pItr = pPriKeyCntsMap->find( rsTblName );
        if( pItr == pPriKeyCntsMap->end( ) )
        {
            // There is no primary key entry for rsTblName so far. Create one ...
            // Use as starting value of the primary key the maximal value so far plus one.
            // DEBUG: std::cout << "In getPtrToPriKeyCnt make new pointer ... " << uiPriKeyCnt << std::endl;
            auto pAutoIncrCntPtr = std::make_shared<std::atomic<PriKeyDefaultType>>( uiPriKeyCnt + 1 );
            pPriKeyCntsMap->emplace( rsTblName, pAutoIncrCntPtr );
            return pAutoIncrCntPtr;
        } // if

        // DEBUG: std::cout << "In getPtrToPriKeyCnt use existing pointer ... " << uiPriKeyCnt << std::endl;
        // There is already an entry for tables of name rsTblName. Deliver the stored pointer to the counter.
        return pItr->second;
    } // method

    /** @brief The function func is executed while protected by the master lock. */
    template <typename F> inline void doSynchronized( F&& func ) const
    {
        std::lock_guard<std::mutex> lock( *pMasterMutex );
        func( );
    } // method
}; // class (SQLDBGlobalSync)

/** @brief Single master object for concurrency synchronization. */
const SQLDBGlobalSync xSQLDBGlobalSync;

/** @brief GLobal literals */
const std::string SCHEMA = "SCHEMA";
const std::string NAME = "NAME";
const std::string DROP_ON_CLOSURE = "DROP_ON_CLOSURE";

/** @brief helper function for the explained query
 * @details some columns return empty strings, in that case we ant to print "NULL".
 */
template <typename Type> inline void printWNull( Type& xEle, fort::char_table& xTable )
{
    xTable << xEle;
} // method

template <> inline void printWNull<std::string>( std::string& xEle, fort::char_table& xTable )
{
    if( xEle.size( ) == 0 )
        xTable << "NULL";
    else
        xTable << xEle;
} // method

template <bool bExplain, typename DBConType> class _QueryExplainer
{
  private:
    const std::string sStmtText; // backup of the query statement text. (used for verbose error reporting)
    const std::string sQueryName;
    // type of the explain query
    using ExplainQType = typename DBConType::template PreparedQuery<int64_t, // id
                                                                    std::string, // select_type
                                                                    std::string, // table
                                                                    std::string, // partitions
                                                                    std::string, // type
                                                                    std::string, // possible_keys
                                                                    std::string, // key
                                                                    std::string, // key_len
                                                                    std::string, // ref
                                                                    uint64_t, // rows
                                                                    double, // filtered
                                                                    std::string // Extra
                                                                    >;
    std::unique_ptr<ExplainQType> pExplainQuery; // pointer to explain statement
    bool bIsExplained; // holds if the query explained itself already

    /** @brief helper for printing the result rows of the explain query */
    template <size_t IDX> struct PrintExplainResults
    {
        template <typename TupleType> bool operator( )( TupleType& xRow, fort::char_table& xTable )
        {
            printWNull( std::get<IDX>( xRow ), xTable ); // print the cell
            return true; // continue the iteration
        } // operator
    }; // struct

    /** @brief print the result rows of the explain query */
    template <typename TupleType> void printExplainRow( TupleType& xRow, fort::char_table& xTable )
    {
        // print cells
        TemplateLoop<std::tuple_size<TupleType>::value, PrintExplainResults>::iterate( xRow, xTable );
        xTable << fort::endr;
    } // method

  public:
    _QueryExplainer( std::shared_ptr<DBConType> pDB, const std::string& rsStmtText,
                     std::string sQueryName = "UnnamedQuery" )
        : sStmtText( rsStmtText ),
          sQueryName( sQueryName ),
          pExplainQuery( bExplain ? std::make_unique<ExplainQType>( pDB, "EXPLAIN " + rsStmtText, json{} ) : nullptr ),
          bIsExplained( false )
    {} // constructor

    template <typename... ArgTypes> void exec( ArgTypes&&... args )
    {
        // explain the query if requested @todo use std::enable_if ?
        if( bExplain /* <- template parameter, so that whole if can be optimized out */ && !bIsExplained )
        {
            bIsExplained = true;
            std::cout << "Explaining " << sQueryName << ":" << std::endl;
            // print query
            std::cout << sStmtText << std::endl;
            // print header
            fort::char_table xTable;
            xTable << fort::header //
                   << "id"
                   << "select_type"
                   << "table"
                   << "partitions"
                   << "type"
                   << "possible_keys"
                   << "key"
                   << "key_len"
                   << "ref"
                   << "rows"
                   << "filtered"
                   << "extra" << fort::endr;
            // Execute statement
            pExplainQuery->bindAndExec( std::forward<ArgTypes>( args )... );
            // Bind outcome of statement execution for later fetching
            pExplainQuery->bindResult( );
            // Fetch all rows and print them
            while( pExplainQuery->fetchNextRow( ) )
                printExplainRow( pExplainQuery->tCellValues, xTable );
            std::cout << xTable.to_string( ) << std::endl;
        } // if
    } // method
}; // class

/** @brief An instance of this class models a single SQL statement. */
template <bool bExplain, typename DBConPtrType> class _SQLStatement
{
  protected:
    DBConPtrType pDB; // Database connection //--
    using DBConType = typename std::remove_reference<decltype( *std::declval<DBConPtrType>( ) )>::type;
    std::unique_ptr<typename DBConType::PreparedStmt> pStmt; // pointer to insert statement
    _QueryExplainer<bExplain, DBConType> xExplainer;

  public:
    _SQLStatement( const _SQLStatement& ) = delete; // no statement copies

    _SQLStatement( DBConPtrType pDB, const std::string& rsStmtText,
                   std::string sStatementName = "UnnamedStatement" ) //--
        : pDB( pDB ),
          pStmt( std::make_unique<typename DBConType::PreparedStmt>( this->pDB, rsStmtText ) ),
          xExplainer( pDB, rsStmtText, sStatementName )
    {} // constructor

    template <typename... ArgTypes> inline void bindDynamic( int uiOffset, ArgTypes&&... args )
    {
        pStmt->bindDynamic( uiOffset, std::forward<ArgTypes>( args )... );
    } // method

    template <int OFFSET, typename... ArgTypes> inline void bindStatic( ArgTypes&&... args )
    {
        pStmt->template bind<OFFSET>( std::forward<ArgTypes>( args )... );
    } // method

    inline DBConPtrType getDBConPtr( )
    {
        return pDB;
    } // method

    /** @brief Connection-safely execute the prepared statement using the DB engine. */
    inline int execWithOutBind( )
    {
        return static_cast<int>( pStmt->exec( ) );
    } // method

    /** @brief Connection-safely execute the prepared statement using the DB engine. */
    template <typename... ArgTypes> inline int exec( ArgTypes&&... args )
    {
        xExplainer.exec( std::forward<ArgTypes>( args )... );
        return static_cast<int>( pStmt->bindAndExec( std::forward<ArgTypes>( args )... ) );
    } // method
}; // class

template <typename DBConType> using SQLStatement = _SQLStatement<false, std::shared_ptr<DBConType>>;
template <typename DBConType> using ExplainedSQLStatement = _SQLStatement<true, std::shared_ptr<DBConType>>;
#define EXPLAINED_QUERY_IN_COMMON


/** @brief An instance of this class models a single SQL query
 *  @details
 *  @param bExplain whether the query shall explain itself on it's first execution
 *  https://stackoverflow.com/questions/24314727/remove-pointer-analog-that-works-for-anything-that-supports-operator/24315093#24315093
 */
#ifdef EXPLAINED_QUERY_IN_COMMON
template <bool bExplain, typename DBConPtrType, typename... ColTypes> class _SQLQueryTmpl
#else
template <typename DBConPtrType, typename... ColTypes> class _SQLQuery
#endif
{
  private:
    DBConPtrType pDB; // Database connection
    // See: https://www.reddit.com/r/cpp_questions/comments/ddgc04/remove_pointert_for_smart_pointers/
    using DBConType = typename std::remove_reference<decltype( *std::declval<DBConPtrType>( ) )>::type;

    std::unique_ptr<typename DBConType::template PreparedQuery<ColTypes...>> pQuery; // pointer to query statement
    const std::string sStmtText; // backup of the query statement text. (used for verbose error reporting)

#ifdef EXPLAINED_QUERY_IN_COMMON
    _QueryExplainer<bExplain, DBConType> xExplainer;
#endif

    template <typename F> void forAllRowsDo( F&& func )
    {
        // Fetch all rows and apply func to all of them
        while( pQuery->fetchNextRow( ) )
            // unpack the tuple and call func (C++17 fold expression)
            STD_APPLY( func, pQuery->tCellValues );
    } // method

  public:
    /** @brief Constructs a query instance using the statement rsStmtText. The query can receive additional parameter by
     *  a JSON  passed via rjConfig.
     */
#ifdef EXPLAINED_QUERY_IN_COMMON
    _SQLQueryTmpl( DBConPtrType pDB, const std::string& rsStmtText, const json& rjConfig = json{},
                   std::string sQueryName = "UnnamedQuery" )
#else
    _SQLQuery( DBConPtrType pDB, const std::string& rsStmtText, const json& rjConfig = json{},
               std::string sQueryName = "UnnamedQuery" )
#endif
        : pDB( pDB ),
          pQuery( std::make_unique<typename DBConType::template PreparedQuery<ColTypes...>>( this->pDB, rsStmtText,
                                                                                             rjConfig ) ),
          sStmtText( rsStmtText )
#ifdef EXPLAINED_QUERY_IN_COMMON
          ,
          xExplainer( pDB, rsStmtText, sQueryName )
#endif
    {} // constructor

    /** @brief Execute the query */
    template <typename... ArgTypes> void exec( ArgTypes&&... args )
    {
#ifdef EXPLAINED_QUERY_IN_COMMON
        // explain the query if requested
        xExplainer.exec( std::forward<ArgTypes>( args )... );
#endif
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
#ifdef POSTGRESQL
        // PostgreSQL requires always reading until EOF or we get an error.
        if( pQuery->fetchNextRow( ) )
            throw std::runtime_error( "SQLQuery: Retrieved more then one row for stmt:" + sStmtText +
                                      "\nArgs: " + dumpArgPack( args... ) );
#endif
        return std::get<COL_NUM>( pQuery->tCellValues );
    } // method

    template <typename... ArgTypes>
    typename std::tuple_element<0, std::tuple<ColTypes...>>::type execAndGetValue( ArgTypes&&... args )
    {
        return execAndGetNthCell<0>( std::forward<ArgTypes>( args )... );
    } // method

    using ColTupleType = std::tuple<ColTypes...>;
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

    /** @brief Delivers a reference to the first cell of the current row. Suitable for queries that deliver a single
     *  value as count(*)-queries etc. */
    inline FstColType& getVal( )
    {
        return std::get<0>( this->get( ) );
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

    /** @brief Delivers the outcome of a query as list of tuples.
     * DEPRECATED!
     */
    template <class... ArgTypes>
    std::vector<typename std::tuple<ColTypes...>> executeAndStoreAllInVector( ArgTypes&&... args )
    {
        // Returned vector that receives column of table
        std::vector<typename std::tuple<ColTypes...>> vRetVec;
        // Execute query, bind, and fetch first row
        this->exec( std::forward<ArgTypes>( args )... );
        // Fetch all rows and add them to the returned vector
        while( pQuery->fetchNextRow( ) )
            vRetVec.push_back( pQuery->tCellValues );
        return vRetVec;
    } // template method


    /** @brief Delivers one column of the outcome of query as vector.
     * DEPRECATED!
     */
    template <size_t COL_NUM, class... ArgTypes>
    std::vector<typename std::tuple_element<COL_NUM, std::tuple<ColTypes...>>::type>
    executeAndStoreInVector( ArgTypes&&... args )
    {
        // Returned vector that receives column of table
        std::vector<typename std::tuple_element<COL_NUM, std::tuple<ColTypes...>>::type> vRetVec;
        // Execute query, bind, and fetch first row
        this->exec( std::forward<ArgTypes>( args )... );
        // Fetch all rows and add them to the returned vector
        while( pQuery->fetchNextRow( ) )
            vRetVec.emplace_back( std::get<COL_NUM>( pQuery->tCellValues ) );
        return vRetVec;
    } // template method

    void explain( )
    {
        // Forward the request to the DB engine.
        pQuery->explain( );
    } // method
}; // class (SQLQuery)

#ifdef EXPLAINED_QUERY_IN_COMMON
/** @brief An instance of this class models a single SQL query. */
template <typename DBCon, typename... ColTypes> using _SQLQuery = _SQLQueryTmpl<false, DBCon, ColTypes...>;

/** @brief An instance of this class models a single SQL query that explains itself on the first usage. */
template <typename DBCon, typename... ColTypes> using _ExplainedSQLQuery = _SQLQueryTmpl<true, DBCon, ColTypes...>;

template <typename DBConType, typename... ColTypes> //
using SQLQuery = _SQLQuery<std::shared_ptr<DBConType>, ColTypes...>;

template <typename DBConType, typename... ColTypes> //
using ExplainedSQLQuery = _ExplainedSQLQuery<std::shared_ptr<DBConType>, ColTypes...>;
#else
/** @brief Single SQL query. */
template <typename DBConType, typename... ColTypes> using SQLQuery = _SQLQuery<std::shared_ptr<DBConType>, ColTypes...>;
#endif


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
const std::string INDEX_TYPE = "INDEX_TYPE";
const std::string INDEX_COLUMNS = "INDEX_COLUMNS";
const std::string WHERE = "WHERE";

// generated columns; see: https://dev.mysql.com/doc/refman/5.7/en/create-table-generated-columns.html
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


/* Tips: https://stackoverflow.com/questions/667643/mysql-batch-updates-in-c
 * FIXME: Change the name to SQLTableView
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

        std::string makeIndexCreateStmt( const std::string& rsTblName, const std::string& rsIdxName,
                                         const std::string& rsNotExistsClause )
        {
            if( jIndexDef.count( INDEX_COLUMNS ) != 1 )
                throw std::runtime_error(
                    "The JSON definition for a SQL index requires one INDEX_COLUMNS item exactly." );
            std::string sCols = jIndexDef[ INDEX_COLUMNS ];

            std::string sIndexType = jIndexDef.count( INDEX_TYPE ) == 0 ? "" : jIndexDef[ INDEX_TYPE ];

            std::string sStmt = std::string( "CREATE " ) +
#ifndef POSTGRESQL // MySQL
                                sIndexType +
#endif
                                " INDEX " + rsNotExistsClause + rsIdxName + " ON " + rsTblName + "(" + sCols + ")";

            // WHERE is not supported by MySQL
            if( jIndexDef.count( WHERE ) )
            {
                if( DBCon::supportsWhereClauseForIndex( ) )
                    sStmt.append( " WHERE " ).append( std::string( jIndexDef[ WHERE ] ) );
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
#ifdef POSTGRESQL
            // PostgreSQL knows the "IF NOT EXISTS" addition.
            pTable->pDB->execSQL( makeIndexCreateStmt( pTable->getTableName( ), sIdxName, "IF NOT EXISTS " ) );
#else /* MySQL */
            // In a pooled environment the creation of indexes should be serialized.
            pTable->pDB->doPoolSafe( [&] {
                // When we check the existence of the index as well as during creation,
                // we require an exclusive lock on the database connection.
                if( !pTable->pDB->indexExists( pTable->getTableName( ), sIdxName ) )
                {
                    pTable->pDB->execSQL( makeIndexCreateStmt( pTable->getTableName( ), sIdxName, "" ) );
                } // if
                else
                {
#ifdef SQL_VERBOSE
                    std::cout << "Index exists already, skip creation." << std::endl;
#endif
                }
            } ); // doPoolSafe
#endif

        } // constructor
    }; // inner class SQL Index

    /** @brief: Inner class for representing view to drop an index.
     * { "INDEX_NAME", "name_of_your_index" } - item is option
     */
    class SQLIndexDropView
    {
        const json jIndexDef; // copy of JSON index definition

        std::string makeIndexDropStmt( const std::string& rsTblName, const std::string& rsIdxName )
        {
            std::string sStmt = std::string( "DROP INDEX " ) +
#ifdef POSTGRESQL
                                "IF NOT EXISTS " +
#endif
                                rsIdxName + " ON " + rsTblName;
            return sStmt;
        } // method

      public:
        SQLIndexDropView( const SQLTable<DBCon, ColTypes...>* pTable,
                          const json& rjIndexDef ) // index definition in JSON)
            : jIndexDef( rjIndexDef ) // copy the index def. to attribute
        {
            // If there is no index name in the JSON, generate one.
            std::string sIdxName = jIndexDef.count( INDEX_NAME ) ? std::string( jIndexDef[ INDEX_NAME ] )
                                                                 : pTable->getTableName( ) + "_idx";

#ifdef POSTGRESQL
            pTable->pDB->execSQL( makeIndexDropStmt( pTable->getTableName( ), sIdxName ) );
#else /* MySQL */
            // In a pooled environment the creation of indexes should be serialized.
            pTable->pDB->doPoolSafe( [&] {
                // When we check the existence of the index as well as during creation,
                // we require an exclusive lock on the database connection.
                if( pTable->pDB->indexExists( pTable->getTableName( ), sIdxName ) )
                {
                    pTable->pDB->execSQL( makeIndexDropStmt( pTable->getTableName( ), sIdxName ) );
                } // if
                else
                {
#ifdef SQL_VERBOSE
                    std::cout << "Index does not exist, skip dropping." << std::endl;
#endif
                } // else
            } ); // doPoolSafe
#endif
        } // constructor
    }; // inner class SQL Index
  public:
    using uiBulkInsertSize = std::integral_constant<size_t, 500>; // number of rows in the insertion buffer

    /** @brief: Implements the concept of bulk inserts for table views.
     *  The bulk inserter always inserts NULL's on columns having type std::nullptr_t.
     */
    template <size_t BUF_SIZE, typename FstType,
              typename... OtherTypes> // InsTypes - either equal to corresponding type in ColTypes or std::nullptr_t
    class SQLBulkInserter
    {
      private:
#ifdef POSTGRESQL
        using InsTupleType =
            typename std::conditional<std::is_same<FstType, PGBigSerial>::value, std::tuple<OtherTypes...>,
                                      std::tuple<FstType, OtherTypes...>>::type;
#else // MySQL
        using InsTupleType = std::tuple<FstType, OtherTypes...>;
#endif
        std::array<InsTupleType, BUF_SIZE> aBuf; // buffer of values inserted by a single op.
        size_t uiInsPos; // index of the next insert into the array
        std::unique_ptr<typename DBCon::PreparedStmt> pBulkInsertStmt; // pointer to prepared bulk insert statement
        std::unique_ptr<typename DBCon::PreparedStmt> pSingleInsertStmt; // required for buffer flush
        std::shared_ptr<std::atomic<PriKeyDefaultType>> pPriKeyCntr; // pointer to auto increment counter

        // If local primary key counting is used, pLocalPriKeyCntr and uiNumRemainPriKey are used for this purpose
        std::unique_ptr<PriKeyDefaultType> pLocalPriKeyCntr; // local primary key counter (nullptr if unused)
        size_t uiNumRemainPriKey; // number of remaining primary keys in batch

#if BULK_INSERT_BY_TUPLE_CATENATION // Nice Haskell like approach, but apparently for most C++ compilers too much  ...
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

#define MAX_COMPILETIME_BIND_N 800

        /** @brief Compile time iteration over the array a for binding each tuple in the array (each row)
         *  @details
         *  This function is only used if the template parameter is N > MAX_COMPILETIME_BIND_N.
         *  In that case we bind the parameters via a runtime loop.
         */
        template <typename T, std::size_t N>
        typename std::enable_if</* condition -> */ ( N > MAX_COMPILETIME_BIND_N ), /* return type -> */ void>::type
        forAllDoBind( const std::array<T, N>& a )
        {
            for( size_t uiOffset = 0; uiOffset < N; uiOffset++ ) // the runtime bind loop
                STD_APPLY( [&]( auto&... args ) //
                           { pBulkInsertStmt->bindDynamic( (int)uiOffset, args... ); },
                           a[ uiOffset ] );
        } // meta

        /** @brief Compile time iteration over the array a for binding each tuple in the array (each row)
         *  @details
         *  This function is only used if the template parameter is N <= MAX_COMPILETIME_BIND_N.
         *  In that case we bind the parameters via a compile-time loop.
         *  @note If this would be compiled with N > 1000 GCC would compile forever and MSVC would throw an exception.
         */
        template <typename T, std::size_t N>
        typename std::enable_if</* condition -> */ ( N <= MAX_COMPILETIME_BIND_N ), /* return type -> */ void>::type
        forAllDoBind( const std::array<T, N>& a )
        {
            forAllDoBindImpl( a, std::make_index_sequence<N>{} );
        } // method


        /** @brief Computes a primary key for a fresh row. */
        inline PriKeyDefaultType getPriKey( )
        {
            assert( pPriKeyCntr != nullptr );

            if( pLocalPriKeyCntr )
            {
                // Keys are assigned from a batch
                if( uiNumRemainPriKey == 0 )
                {
                    // All locals keys are exhausted, get a new batch of BUF_SIZE many keys
                    uiNumRemainPriKey = BUF_SIZE;
                    *pLocalPriKeyCntr = pPriKeyCntr->fetch_add( BUF_SIZE );
                } // if
                // Return a key from the current batch
                uiNumRemainPriKey--;
                return ( *pLocalPriKeyCntr )++;
            } // if
            else
                // Each key is 'globally' requested.
                return pPriKeyCntr->fetch_add( 1 );
        } // method

        /** @brief If the buffer is full, insert the buffer content into the table by executing a bulk insert
         * statement. Design alternative: iterate over the array and bind each tuple separately. Discussion: The
         * latter approach requires more code but avoids the template instantiation problems.
         */
        inline void bulkInsertIfBufFull( )
        {
            // If the buffer is full, do a bulk insert
            if( ++uiInsPos >= BUF_SIZE )
            {
                // The DB operations does not inspects uiInsPos anymore. However, they can throw an exceptions.
                // In order to avoid the crash scenario described in flush(), uiInsPos must be set zero before
                // doing any DB ops.
                uiInsPos = 0; // reset counter
                // Write the buffer to the database
                // Step 1: Iterates over the array during compile time and binds them.
                forAllDoBind( aBuf );

                // Step 2: After binding all arguments, actually execute the insert statement.
                pBulkInsertStmt->exec( );
            } // if
        } // method

        /** @brief Store a single row in the buffer and return the primary key. */
        inline PriKeyDefaultType storeInBuf( Identity<PriKeyDefaultType> _, const OtherTypes&... args )
        {
            // If required book a new batch of primary keys.
            PriKeyDefaultType uiPriKey = getPriKey( );
            aBuf[ uiInsPos ] = InsTupleType( uiPriKey, args... );
            return uiPriKey;
        } // method

        /** @brief Store a single row in the buffer and return a nullptr. */
        inline std::nullptr_t storeInBuf( Identity<std::nullptr_t> _, const OtherTypes&... args )
        {
            aBuf[ uiInsPos ] = InsTupleType( nullptr, args... );
            return nullptr;
        } // method

#ifdef POSTGRESQL
        /** @brief Store a single row in the buffer and return a nullptr. */
        inline PGBigSerial storeInBuf( Identity<PGBigSerial> _, const OtherTypes&... args )
        {
            aBuf[ uiInsPos ] = InsTupleType( args... );
            return PGBigSerial( 0 );
        } // method
#endif // POSTGRESQL
#endif // alternative approach

      public:
        /** @brief Construct bulk-inserter.
         *  The table must be part of the database already, or the construction fails.
         *  @details
         *  If bReserveBatchOfPriKeys is true, the bulk inserter reserves a batch of keys, where the size of the batch
         *  is equal to the size of the buffer. The bulk inserter uses the keys in this batch for assigning keys in
         *  strictly increasing order. This helps the database regarding the insertion of keys in the B-tree (index) and
         *  speeds up the bulk inserter in a multi-threaded environment.
         *  Implementation detail: Because we can not guarantee the life of host table for all of the life of the bulk
         *  inserter, we should not keep a reference to the table. Instead we copy the required info of the table.
         */
        SQLBulkInserter( const SQLTable& rxHost, // host table
                         std::shared_ptr<std::atomic<PriKeyDefaultType>> pPriKeyCntr = nullptr,
                         bool bReserveBatchOfPriKeys = false )
            : uiInsPos( 0 ),
              pBulkInsertStmt( // compile bulk insert stmt
                  std::make_unique<typename DBCon::PreparedStmt>( rxHost.pDB, rxHost.makeInsertStmt( BUF_SIZE ) ) ),
              pSingleInsertStmt( // compile single insert stmt
                  std::make_unique<typename DBCon::PreparedStmt>( rxHost.pDB, rxHost.makeInsertStmt( 1 ) ) ),
              pPriKeyCntr( pPriKeyCntr ),
              pLocalPriKeyCntr( bReserveBatchOfPriKeys ? std::make_unique<PriKeyDefaultType>( 0 ) : nullptr ),
              uiNumRemainPriKey( 0 ) // initially, the batch is empty
        {
            static_assert( sizeof...( OtherTypes ) + 1 == sizeof...( ColTypes ) );
            pBulkInsertStmt->setArgsMultiplicator( BUF_SIZE );
        } // constructor

        /**
         * @brief get a bulk inserter directly from a connection
         * @details
         * This constructs the respective table and then immediately destructs it again.
         * That way we make sure the table exists in the DB.
         */
        SQLBulkInserter( std::shared_ptr<DBCon> pConnection );

        /** @brief Inserts a row into the table via a bulk-insert approach.
         *  Reasonable addition: insert using moves.
         */
        inline void insert( const FstType& fstArg, const OtherTypes&... otherArgs )
        {
            assert( pPriKeyCntr == nullptr );

            aBuf[ uiInsPos ] = InsTupleType( fstArg, otherArgs... ); // This triggers a copy...
            bulkInsertIfBufFull( ); // write buffer to DB if full
        } // method

        /** @brief Inserts a row into the table via a bulk-insert approach for tables with automatic primary key. */
        inline FstType insert( const OtherTypes&... args )
        {
            // First column is expected to keep a primary key. In the case of AUTO_INCREMENT, FstType is
            // std::nullptr_t. In the case of a primary key that is incremented by the library, FstType is
            // PriKeyDefaultType.
            static_assert( ( std::is_same<FstType, std::nullptr_t>::value ) ||
#ifdef POSTGRESQL
                               ( std::is_same<FstType, PGBigSerial>::value ) ||
#endif
                               ( std::is_same<FstType, PriKeyDefaultType>::value ),
                           "FstType must be either std::nullptr_t or PriKeyDefaultType" );
            // pAutoIncrCounter must be initialized if FstType is PriKeyDefaultType
            assert( !( std::is_same<FstType, PriKeyDefaultType>::value && ( pPriKeyCntr == nullptr ) ) );

            FstType uiRetVal = storeInBuf( Identity<FstType>( ), args... );
            bulkInsertIfBufFull( ); // write buffer to DB if full
            return uiRetVal;
        } // method

        /** @brief Flush a non-full buffer to the database. */
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

        /** @brief Destructor flushes the buffer. */
        ~SQLBulkInserter( )
        {
            // Throwing an exception in a destructor results in undefined behavior.
            // Therefore, we swallow these exceptions and report them via std:cerr.
            doNoExcept( [this] { this->flush( ); }, "Exception in ~SQLBulkInserter:" );
        } // destructor
    }; // class (SQLBulkInserter)
#if 0
#define NUM_CON 1
    /** @brief: Implements the concept of bulk inserts for table views.
     *  The bulk inserter always inserts NULL's on columns having type std::nullptr_t.
     */
    template <typename SQLDBConPoolType, size_t BUF_SIZE, typename FstType,
              typename... OtherTypes> // InsTypes - either equal to corresponding type in ColTypes or std::nullptr_t
    class AsyncSQLBulkInserter
    {
      private:
        using InsTupleType = std::tuple<FstType, OtherTypes...>;

        std::shared_ptr<SQLDBConPoolType> pxConPool; // Connection pool the delivers two dedicated connections.
        // std::array<InsTupleType, BUF_SIZE> aBuf; // buffer of values inserted by a single op.
        std::array<std::array<InsTupleType, BUF_SIZE>, NUM_CON> aBuf;
        size_t uiInsPos; // index of the next insert into the array

        std::array<std::unique_ptr<typename DBCon::PreparedStmt>, NUM_CON>
            pBulkInsertStmt; // pointer to prepared bulk insert statement
        std::array<std::unique_ptr<typename DBCon::PreparedStmt>, NUM_CON>
            pSingleInsertStmt; // required for buffer flush

        std::shared_ptr<std::atomic<PriKeyDefaultType>> pPriKeyCntr; // pointer to auto increment counter

        // If local primary key counting is used, pLocalPriKeyCntr and uiNumRemainPriKey are used for this purpose
        std::unique_ptr<PriKeyDefaultType> pLocalPriKeyCntr; // local primary key counter (nullptr if unused)
        size_t uiNumRemainPriKey; // number of remaining primary keys in batch

        /** @brief Unpack the tuple to an argument pack and bind using this pack.
         *  (Part of bulkInsert)
         */
        template <int OFFSET, typename... Args> void doSingleBind( const std::tuple<Args...>& rTpl )
        {
            STD_APPLY( [ & ]( auto&... args ) //
                       { pBulkInsertStmt[ 0 ]->template bind<OFFSET>( args... ); },
                       rTpl );
        } // method

        /** @brief Implementation part of forAllDoBind */
        template <typename T, std::size_t N, std::size_t... Idx>
        void forAllDoBindImpl( const std::array<T, N>& a, std::index_sequence<Idx...> )
        {
            // Requires C++17 folded expression support ...
            ( doSingleBind<Idx>( a[ Idx ] ), ... );
        } // meta

        //-- #define MAX_COMPILETIME_BIND_N 800

        /** @brief Compile time iteration over the array a for binding each tuple in the array (each row)
         *  @details
         *  This function is only used if the template parameter is N > MAX_COMPILETIME_BIND_N.
         *  In that case we bind the parameters via a runtime loop.
         */
        //-- template <typename T, std::size_t N>
        //-- typename std::enable_if</* condition -> */ ( N > MAX_COMPILETIME_BIND_N ), /* return type -> */ void>::type
        //-- forAllDoBind( const std::array<T, N>& a )
        //-- {
        //--     for( size_t uiOffset = 0; uiOffset < N; uiOffset++ ) // the runtime bind loop
        //--         STD_APPLY( [ & ]( auto&... args ) //
        //--                    { pBulkInsertStmt[0]->bindDynamic( (int)uiOffset, args... ); },
        //--                    a[ uiOffset ] );
        //-- } // meta

        /** @brief Compile time iteration over the array a for binding each tuple in the array (each row)
         *  @details
         *  This function is only used if the template parameter is N <= MAX_COMPILETIME_BIND_N.
         *  In that case we bind the parameters via a compile-time loop.
         *  @note If this would be compiled with N > 1000 GCC would compile forever and MSVC would throw an exception.
         */
        template <typename T, std::size_t N>
        typename std::enable_if</* condition -> */ ( N <= MAX_COMPILETIME_BIND_N ), /* return type -> */ void>::type
        forAllDoBind( const std::array<T, N>& a )
        {
            forAllDoBindImpl( a, std::make_index_sequence<N>{ } );
        } // method


        /** @brief Computes a primary key for a fresh row. */
        inline PriKeyDefaultType getPriKey( )
        {
            assert( pPriKeyCntr != nullptr );

            if( pLocalPriKeyCntr )
            {
                // Keys are assigned from a batch
                if( uiNumRemainPriKey == 0 )
                {
                    // All locals keys are exhausted, get a new batch of BUF_SIZE many keys
                    uiNumRemainPriKey = BUF_SIZE;
                    *pLocalPriKeyCntr = pPriKeyCntr->fetch_add( BUF_SIZE );
                } // if
                // Return a key from the current batch
                uiNumRemainPriKey--;
                return ( *pLocalPriKeyCntr )++;
            } // if
            else
                // Each key is 'globally' requested.
                return pPriKeyCntr->fetch_add( 1 );
        } // method

        std::array<std::future<void>, NUM_CON> aFutures; // array keeping futures for bulk inserts
        std::array<int, NUM_CON> aThreadId; // array keeping ids of dedicated connections
        std::string sInsStmtTxt;

        /** @brief If the buffer is full, insert the buffer content into the table by executing a bulk insert
         * statement. Design alternative: iterate over the array and bind each tuple separately. Discussion: The
         * latter approach requires more code but avoids the template instantiation problems.
         */
        inline void bulkInsertIfBufFull( )
        {
            // If the buffer is full, do a bulk insert
            if( ++uiInsPos >= BUF_SIZE )
            {
#if 0
                // The DB operations does not inspects uiInsPos anymore. However, they can throw an exceptions.
                // In order to avoid the crash scenario described in flush(), uiInsPos must be set zero before
                // doing any DB ops.
                uiInsPos = 0; // reset counter
                // Write the buffer to the database
                // Step 1: Iterates over the array during compile time and binds them.
                forAllDoBind( aBuf[ 0 ] );

                // Step 2: After binding all arguments, actually execute the insert statement.
                pBulkInsertStmt[ 0 ]->exec( );
#endif
                // Guarantee that the previous exec for the selected connection finished already.
                if( aFutures[ 0 ].valid( ) )
                    aFutures[ 0 ].get( );

                // The DB operations does not inspects uiInsPos anymore. However, they can throw an exceptions.
                // In order to avoid the crash scenario described in flush(), uiInsPos must be set zero before
                // doing any DB ops.
                uiInsPos = 0; // reset counter
                // Write the buffer to the database
                // Step 1: Iterates over the array during compile time and binds them.
                // Run and wait until the finished (until future is available).
                auto xFut = pxConPool->enqueue(
                    aThreadId[ 0 ],
                    [ & ]( auto pCon, int id ) {
                        pBulkInsertStmt[0] =
                            std::make_unique<typename DBCon::PreparedStmt>(pCon, sInsStmtTxt);
                        forAllDoBind( aBuf[ 0 ] );
                        pBulkInsertStmt[ 0 ]->exec( );
                    },
                    0 );
                xFut.get( );

                // Step 2: After binding all arguments, actually execute the insert statement.
                // aFutures[ 0 ] = pxConPool->enqueue( aThreadId[ 0 ],
                //                                    [ & ]( auto pCon, int id ) { pBulkInsertStmt[ 0 ]->exec( ); }, 0
                //                                    );
                // aFutures[ 0 ].get( );
            } // if
        } // method

        /** @brief Store a single row in the buffer and return the primary key. */
        inline PriKeyDefaultType storeInBuf( Identity<PriKeyDefaultType> _, const OtherTypes&... args )
        {
            // If required book a new batch of primary keys.
            PriKeyDefaultType uiPriKey = getPriKey( );
            aBuf[ 0 ][ uiInsPos ] = InsTupleType( uiPriKey, args... );
            return uiPriKey;
        } // method

        /** @brief Store a single row in the buffer and return a nullptr. */
        inline std::nullptr_t storeInBuf( Identity<std::nullptr_t> _, const OtherTypes&... args )
        {
            aBuf[ 0 ][ uiInsPos ] = InsTupleType( nullptr, args... );
            return nullptr;
        } // method

      public:
        /** @brief Construct bulk-inserter.
         *  The table must be part of the database already, or the construction fails.
         *  @details
         *  If bReserveBatchOfPriKeys is true, the bulk inserter reserves a batch of keys, where the size of the batch
         *  is equal to the size of the buffer. The bulk inserter uses the keys in this batch for assigning keys in
         *  strictly increasing order. This helps the database regarding the insertion of keys in the B-tree (index) and
         *  speeds up the bulk inserter in a multi-threaded environment.
         *  Implementation detail: Because we can not guarantee the life of host table for all of the life of the bulk
         *  inserter, we should not keep a reference to the table. Instead we copy the required info of the table.
         */
        AsyncSQLBulkInserter( std::shared_ptr<SQLDBConPoolType> pxConPool, // connection pool
                              const SQLTable& rxHost, // host table
                              std::shared_ptr<std::atomic<PriKeyDefaultType>> pPriKeyCntr = nullptr,
                              bool bReserveBatchOfPriKeys = false )
            : pxConPool( pxConPool ), // store local copy of pointer to con pool
              uiInsPos( 0 ),
              pSingleInsertStmt( // compile single insert stmt
                  { std::make_unique<typename DBCon::PreparedStmt>( rxHost.pDB, rxHost.makeInsertStmt( 1 ) ) } ),
              pPriKeyCntr( pPriKeyCntr ),
              pLocalPriKeyCntr( bReserveBatchOfPriKeys ? std::make_unique<PriKeyDefaultType>( 0 ) : nullptr ),
              uiNumRemainPriKey( 0 ), // initially, the batch is empty
              aFutures( { std::future<void>( ) } ),
              sInsStmtTxt( rxHost.makeInsertStmt( BUF_SIZE ) )
        {
            static_assert( sizeof...( OtherTypes ) + 1 == sizeof...( ColTypes ) );

            aThreadId[ 0 ] = pxConPool->getDedicatedConId( );


            auto xFut = pxConPool->enqueue( //
                aThreadId[ 0 ],
                [ & ]( auto pCon, int id ) {
                    pBulkInsertStmt[ 0 ] =
                        std::make_unique<typename DBCon::PreparedStmt>( pCon, sInsStmtTxt );
                },
                0 );
            xFut.get( );
        } // constructor

        /**
         * @brief get a bulk inserter directly from a connection
         * @details
         * This constructs the respective table and then immediately destructs it again.
         * That way we make sure the table exists in the DB.
         */
        // AsyncSQLBulkInserter(std::shared_ptr<DBCon> pConnection);

        /** @brief Inserts a row into the table via a bulk-insert approach.
         *  Reasonable addition: insert using moves.
         */
        inline void insert( const FstType& fstArg, const OtherTypes&... otherArgs )
        {
            assert( pPriKeyCntr == nullptr );

            aBuf[ 0 ][ uiInsPos ] = InsTupleType( fstArg, otherArgs... ); // This triggers a copy...
            bulkInsertIfBufFull( ); // write buffer to DB if full
        } // method

        /** @brief Inserts a row into the table via a bulk-insert approach for tables with automatic primary key. */
        inline FstType insert( const OtherTypes&... args )
        {
            // First column is expected to keep a primary key. In the case of AUTO_INCREMENT, FstType is
            // std::nullptr_t. In the case of a primary key that is incremented by the library, FstType is
            // PriKeyDefaultType.
            static_assert( ( std::is_same<FstType, std::nullptr_t>::value ) ||
                               ( std::is_same<FstType, PriKeyDefaultType>::value ),
                           "FstType must be either std::nullptr_t or PriKeyDefaultType" );
            // pAutoIncrCounter must be initialized if FstType is PriKeyDefaultType
            assert( !( std::is_same<FstType, PriKeyDefaultType>::value && ( pPriKeyCntr == nullptr ) ) );

            FstType uiRetVal = storeInBuf( Identity<FstType>( ), args... );
            bulkInsertIfBufFull( ); // write buffer to DB if full
            return uiRetVal;
        } // method

        /** @brief Flush a non-full buffer to the database. */
        inline void flush( )
        {
            // uiInsPos must be set zero before calling bindAndExec, because bindAndExec could throw an exception
            // which in turn triggers a destructor call. However, in the case of an exception, the flush in the
            // destructor should not lead to an additional call of bindAndExec, because it would lead to a crash.
            auto uiInsPosBackup = uiInsPos;
            uiInsPos = 0; // reset insert position counter
            for( size_t uiCount = 0; uiCount < uiInsPosBackup; uiCount++ )
                // Write the current tuple in the array to the DB
                STD_APPLY( [ & ]( auto&... args ) { pSingleInsertStmt[ 0 ]->bindAndExec( args... ); },
                           aBuf[ 0 ][ uiCount ] );
        } // method

        /** @brief Destructor flushes the buffer. */
        ~AsyncSQLBulkInserter( )
        {
            // Throwing an exception in a destructor results in undefined behavior.
            // Therefore, we swallow these exceptions and report them via std:cerr.
            doNoExcept( [ this ] { this->flush( ); }, "Exception in ~SQLBulkInserter:" );
        } // destructor
    }; // class (AsyncSQLBulkInserter)
#endif
#if 1
    /** @brief: Implements the concept of bulk inserts for table views.
     *  The bulk inserter always inserts NULL's on columns having type std::nullptr_t.
     */
#define NUM_CON 2
    template <typename SQLDBConPoolType, size_t BUF_SIZE, typename FstType,
              typename... OtherTypes> // InsTypes - either equal to corresponding type in ColTypes or std::nullptr_t
    class AsyncSQLBulkInserter
    {
      private:
        using InsTupleType = std::tuple<FstType, OtherTypes...>;
        using ConIdType = int; // TODO: get this type from the pool

        std::shared_ptr<SQLDBConPoolType> pxConPool; // Connection pool the delivers two dedicated connections.
        std::array<std::array<InsTupleType, BUF_SIZE>, NUM_CON> aBufPair; // buffer of values inserted by a single op.
        InsTupleType* pCurrBuf; // pointer to active buffer in buffer pair
        size_t uiSelConId; // index of currently selected dedicated connection.
        size_t uiInsPos; // index of the next insert into the array
        std::array<std::unique_ptr<SQLStatement<DBCon>>, NUM_CON>
            apBulkInsertStmt; // pointer to prepared bulk insert statement
        std::array<std::unique_ptr<SQLStatement<DBCon>>, NUM_CON> apSingleInsertStmt; // required for buffer flush
        std::shared_ptr<std::atomic<PriKeyDefaultType>> pPriKeyCntr; // pointer to auto increment counter

        // If local primary key counting is used, pLocalPriKeyCntr and uiNumRemainPriKey are used for this purpose
        std::unique_ptr<PriKeyDefaultType> pLocalPriKeyCntr; // local primary key counter (nullptr if unused)
        size_t uiNumRemainPriKey; // number of remaining primary keys in batch

        /** @brief Unpack the tuple to an argument pack and bind using this pack.
         *  (Part of bulkInsert)
         */
        template <int OFFSET, typename... Args> void doSingleBind( size_t iConId, const std::tuple<Args...>& rTpl )
        {
            STD_APPLY( [&]( auto&... args ) //
                       { apBulkInsertStmt[ iConId ]->template bindStatic<OFFSET>( args... ); },
                       rTpl );
        } // method

        /** @brief Implementation part of forAllDoBind */
        template <typename T, std::size_t N, std::size_t... Idx>
        void forAllDoBindImpl( size_t iConId, const std::array<T, N>& a, std::index_sequence<Idx...> )
        {
            // Requires C++17 folded expression support ...
            ( doSingleBind<Idx>( iConId, a[ Idx ] ), ... );
        } // meta

        //--template <typename T, std::size_t... Idx> void forAllDoBindImpl( int id, T* a, std::index_sequence<Idx...> )
        //--{
        //--    // Requires C++17 folded expression support ...
        //--    ( doSingleBind<Idx>( id, a[ Idx ] ), ... );
        //--} // meta

        /** @brief Compile time iteration over the array a for binding each tuple in the array (each row)
         *  @details
         *  This function is only used if the template parameter is N > MAX_COMPILETIME_BIND_N.
         *  In that case we bind the parameters via a runtime loop.
         */
        template <typename T, std::size_t N>
        typename std::enable_if</* condition -> */ ( N > MAX_COMPILETIME_BIND_N ), /* return type -> */ void>::type
        forAllDoBind( int id, const std::array<T, N>& a )
        {
            for( size_t uiOffset = 0; uiOffset < N; uiOffset++ ) // the runtime bind loop
                STD_APPLY( [&]( auto&... args ) //
                           { apBulkInsertStmt[ id ]->bindDynamic( (int)uiOffset, args... ); },
                           a[ uiOffset ] );
        } // meta

        //-- template <typename T, std::size_t N>
        //-- typename std::enable_if</* condition -> */ ( N > MAX_COMPILETIME_BIND_N ), /* return type -> */ void>::type
        //-- forAllDoBind( int id, T* a )
        //-- {
        //--     for( size_t uiOffset = 0; uiOffset < N; uiOffset++ ) // the runtime bind loop
        //--         STD_APPLY( [ & ]( auto&... args ) //
        //--                    { apBulkInsertStmt[ id ]->bindDynamic( (int)uiOffset, args... ); },
        //--                    a[ uiOffset ] );
        //-- } // meta

        /** @brief Compile time iteration over the array a for binding each tuple in the array (each row)
         *  @details
         *  This function is only used if the template parameter is N <= MAX_COMPILETIME_BIND_N.
         *  In that case we bind the parameters via a compile-time loop.
         *  @note If this would be compiled with N > 1000 GCC would compile forever and MSVC would throw an exception.
         */
        template <typename T, std::size_t N>
        typename std::enable_if</* condition -> */ ( N <= MAX_COMPILETIME_BIND_N ), /* return type -> */ void>::type
        forAllDoBind( size_t iConId, const std::array<T, N>& a )
        {
            forAllDoBindImpl( iConId, a, std::make_index_sequence<N>{} );
        } // method

        //--template <typename T, std::size_t N>
        //--typename std::enable_if</* condition -> */ ( N <= MAX_COMPILETIME_BIND_N ), /* return type -> */ void>::type
        //--forAllDoBind( int id, T* a )
        //--{
        //--    forAllDoBindImpl( id, a, std::make_index_sequence<N>{ } );
        //--} // method


        /** @brief Computes a primary key for a fresh row. */
        inline PriKeyDefaultType getPriKey( )
        {
            assert( pPriKeyCntr != nullptr );

            if( pLocalPriKeyCntr )
            {
                // Keys are assigned from a batch
                if( uiNumRemainPriKey == 0 )
                {
                    // All locals keys are exhausted, get a new batch of BUF_SIZE many keys
                    uiNumRemainPriKey = BUF_SIZE;
                    *pLocalPriKeyCntr = pPriKeyCntr->fetch_add( BUF_SIZE );
                } // if
                // Return a key from the current batch
                uiNumRemainPriKey--;
                return ( *pLocalPriKeyCntr )++;
            } // if
            else
                // Each key is 'globally' requested.
                return pPriKeyCntr->fetch_add( 1 );
        } // method


        std::array<std::future<void>, NUM_CON> aFutures; // array keeping futures for bulk inserts
        std::array<int, NUM_CON> aThreadIds; // array keeping ids of dedicated connections

        /** @brief If the buffer is full, insert the buffer content into the table by executing a bulk insert
         * statement. Design alternative: iterate over the array and bind each tuple separately. Discussion: The
         * latter approach requires more code but avoids the template instantiation problems.
         */
        inline void bulkInsertIfBufFull( )
        {
            // If the buffer is full, do a bulk insert
            if( ++uiInsPos >= BUF_SIZE )
            {
                // The DB operations does not inspects uiInsPos anymore. However, they can throw an exceptions.
                // In order to avoid the crash scenario described in flush(), uiInsPos must be set zero before
                // doing any DB ops.
                uiInsPos = 0; // reset counter
                // Write the buffer to the database
                // Step 1: Iterates over the array during compile time and binds them.
                // Run and wait until the finished (until future is available).
                pxConPool->run(
                    aThreadIds[ uiSelConId ],
                    [&]( auto pCon, size_t iConId ) { this->forAllDoBind( iConId, aBufPair[ iConId ] ); },
                    uiSelConId );

                // Step 2: After binding all arguments, actually execute the insert statement.
                aFutures[ uiSelConId ] = pxConPool->enqueue(
                    aThreadIds[ uiSelConId ],
                    [&]( auto pCon, size_t iConId ) { apBulkInsertStmt[ iConId ]->execWithOutBind( ); },
                    uiSelConId );
                aFutures[ uiSelConId ].get( );
                // Step 3: Switch to the opposite connection.
                uiSelConId = !uiSelConId; // alters between zero and one
                pCurrBuf = aBufPair[ uiSelConId ].data( ); // set the buffer pointer to the opposite connection
            } // if
        } // method

        /** @brief Wait until the current connection finished its insertion into the DB */
        void waitUntilInsDone( )
        {
            // if( aFutures[ uiSelConId ].valid( ) )
            //     aFutures[ uiSelConId ].get( );
        } // method

        /** @brief Store a single row in the buffer and return the primary key. */
        inline PriKeyDefaultType storeInBuf( Identity<PriKeyDefaultType> _, const OtherTypes&... args )
        {
            assert( pCurrBuf == aBufPair[ uiSelConId ].data( ) );
            // Guarantee that the previous exec for the selected connection finished already.
            // If required book a new batch of primary keys.
            PriKeyDefaultType uiPriKey = getPriKey( );
            pCurrBuf[ uiInsPos ] = InsTupleType( uiPriKey, args... );
            return uiPriKey;
        } // method

        /** @brief Store a single row in the buffer and return a nullptr. */
        inline std::nullptr_t storeInBuf( Identity<std::nullptr_t> _, const OtherTypes&... args )
        {
            // Guarantee that the previous exec for the selected connection finished already.
            pCurrBuf[ uiInsPos ] = InsTupleType( nullptr, args... );
            return nullptr;
        } // method

      public:
        /** @brief Construct bulk-inserter.
         *  The table must be part of the database already, or the construction fails.
         *  @details
         *  If bReserveBatchOfPriKeys is true, the bulk inserter reserves a batch of keys, where the size of the
         * batch is equal to the size of the buffer. The bulk inserter uses the keys in this batch for assigning
         * keys in strictly increasing order. This helps the database regarding the insertion of keys in the B-tree
         * (index) and speeds up the bulk inserter in a multi-threaded environment. Implementation detail: Because
         * we can not guarantee the life of host table for all of the life of the bulk inserter, we should not keep
         * a reference to the table. Instead we copy the required info of the table.
         */
        AsyncSQLBulkInserter( std::shared_ptr<SQLDBConPoolType> pxConPool,
                              const SQLTable& rxHost, // host table
                              std::shared_ptr<std::atomic<PriKeyDefaultType>> pPriKeyCntr = nullptr,
                              bool bReserveBatchOfPriKeys = false )
            : pxConPool( pxConPool ), // store local copy of pointer to con pool
              pCurrBuf( aBufPair[ 0 ].data( ) ),
              uiSelConId( 0 ), // todo: move before pCurrBuf
              uiInsPos( 0 ),
              pPriKeyCntr( pPriKeyCntr ),
              pLocalPriKeyCntr( bReserveBatchOfPriKeys ? std::make_unique<PriKeyDefaultType>( 0 ) : nullptr ),
              uiNumRemainPriKey( 0 ), // initially, the batch is empty
              aFutures( {std::future<void>( ), std::future<void>( )} ) // initialzied futures

        {
            static_assert( sizeof...( OtherTypes ) + 1 == sizeof...( ColTypes ) );

            // Get the ids for two dedicated connections that will be used for bulk inserts
            aThreadIds[ 0 ] = pxConPool->getDedicatedConId( );
            aThreadIds[ 1 ] = pxConPool->getDedicatedConId( );

            // Create the required insertion statements within both dedicated connections for later use.
            for( size_t uiItr = 0; uiItr < NUM_CON; uiItr++ )
            {
                pxConPool->run( (int)aThreadIds[ uiItr ],
                                [&]( auto pCon, size_t id ) {
                                    // Compile bulk insert stmt
                                    apBulkInsertStmt[ id ] = std::make_unique<SQLStatement<DBCon>>(
                                        pCon, rxHost.makeInsertStmt( BUF_SIZE ) );
                                    // compile single insert stmt
                                    apSingleInsertStmt[ id ] =
                                        std::make_unique<SQLStatement<DBCon>>( pCon, rxHost.makeInsertStmt( 1 ) );
                                },
                                uiItr );
            } // for
        } // constructor

        /**
         * @brief get a bulk inserter directly from a connection
         * @details
         * This constructs the respective table and then immediately destructs it again.
         * That way we make sure the table exists in the DB.
         */
        // SQLBulkInserter( std::shared_ptr<DBCon> pConnection );

        /** @brief Inserts a row into the table via a bulk-insert approach.
         *  Reasonable addition: insert using moves.
         */
        inline void insert( const FstType& fstArg, const OtherTypes&... otherArgs )
        {
            assert( pPriKeyCntr == nullptr );

            waitUntilInsDone( ); // wait until the current connection is available
            pCurrBuf[ uiInsPos ] = InsTupleType( fstArg, otherArgs... ); // This triggers a copy...
            bulkInsertIfBufFull( ); // write buffer to DB if full
        } // method

        /** @brief Inserts a row into the table via a bulk-insert approach for tables with
         *  automatic primary key.
         */
        inline FstType insert( const OtherTypes&... args )
        {
            // First column is expected to keep a primary key. In the case of AUTO_INCREMENT, FstType is
            // std::nullptr_t. In the case of a primary key that is incremented by the library, FstType is
            // PriKeyDefaultType.
            static_assert( ( std::is_same<FstType, std::nullptr_t>::value ) ||
                               ( std::is_same<FstType, PriKeyDefaultType>::value ),
                           "FstType must be either std::nullptr_t or PriKeyDefaultType" );
            // pAutoIncrCounter must be initialized if FstType is PriKeyDefaultType
            assert( !( std::is_same<FstType, PriKeyDefaultType>::value && ( pPriKeyCntr == nullptr ) ) );

            waitUntilInsDone( ); // wait until the current connection is available
            FstType uiRetVal = storeInBuf( Identity<FstType>( ), args... );
            bulkInsertIfBufFull( ); // write buffer to DB if full
            return uiRetVal;
        } // method

        /** @brief Flush a non-full buffer to the database. */
        inline void flush( )
        {
            // Guarantee that both connections finished their exec.
            for( size_t uiItr = 0; uiItr < NUM_CON; uiItr++ )
            {
                if( aFutures[ uiItr ].valid( ) )
                    aFutures[ uiItr ].get( );
            } // for

            // uiInsPos must be set zero before calling bindAndExec, because bindAndExec could throw an
            // exception which in turn triggers a destructor call. However, in the case of an exception, the
            // flush in the destructor should not lead to an additional call of bindAndExec, because it
            // would lead to a crash.
            auto uiInsPosBackup = uiInsPos;
            uiInsPos = 0; // reset insert position counter

            // Write the remaining inserts using connection 0.
            pxConPool->run( (int)aThreadIds[ 0 ], [&]( auto pCon ) {
                for( size_t uiCount = 0; uiCount < uiInsPosBackup; uiCount++ )
                    // Write the current tuple in the array to the DB
                    STD_APPLY( [&]( auto&... args ) { apSingleInsertStmt[ 0 ]->exec( args... ); },
                               pCurrBuf[ uiCount ] );
            } );
        } // method

        /** @brief Destructor flushes the buffer. */
        ~AsyncSQLBulkInserter( )
        {
            // Throwing an exception in a destructor results in undefined behavior.
            // Therefore, we swallow these exceptions and report them via std:cerr.
            doNoExcept( [this] { this->flush( ); }, "Exception in ~SQLBulkInserter:" );
        } // destructor
    }; // class (AsyncSQLBulkInserter)
#endif
#if 0
    /** @brief Like a standard SQLBulkInserter, but for the streaming a CSV-file is used in between. This
     * approach is sometimes faster than the standard SQLBulkInserter. IMPORTANT NOTICE: Call release() before a
     * SQLFileBulkInserter runs out of scope. In this case you get an exception, if something goes wrong. If you
     * let do the release-job automatically by the destructor all exception are swallowed and you won't get any
     * feedback about problems.
     */
    template <size_t BUF_SIZE,
              typename... InsTypes> // InsTypes - either equal to corresponding type in ColTypes or
                                    // std::nullptr_t
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
            aBuf[ uiInsPos ] = std::tuple<InsTypes...>( args... );
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
#endif
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
        // SQL code for regular columns
        for( size_t iItr = 0; iItr < this->rjTableCols.size( ); )
        {
            auto& rjCol = this->rjTableCols[ iItr ]; // current column in jTableDef
            sStmt
                .append( std::string( rjCol[ COLUMN_NAME ] ) ) // column name
                .append( " " )
                .append( vSQLColumnTypes[ iItr ] ); // column type
            if( rjCol.count( CONSTRAINTS ) )
                // CONSTRAINTS must be plain text describing all constraints
                sStmt.append( " " ).append( std::string( rjCol[ CONSTRAINTS ] ) );
            // insert separating comma
            if( ++iItr < this->rjTableCols.size( ) )
                sStmt.append( ", " );
        } // for

        // SQL code for generated columns if there are any defined
        if( this->jTableDef.count( GENERATED_COLUMNS ) )
            for( size_t iItr = 0; iItr < this->jTableDef[ GENERATED_COLUMNS ].size( ); iItr++ )
            {
                // get the json column entry
                auto& rjCol = this->jTableDef[ GENERATED_COLUMNS ][ iItr ]; // current column in jTableDef

                // check that each column is complete
                if( rjCol.count( COLUMN_NAME ) == 0 )
                    throw std::runtime_error(
                        std::string( "COLUMN_NAME is missing for a GENERATED_COLUMNS in the table " )
                            .append( getTableName( ) ) );
                if( rjCol.count( TYPE ) == 0 )
                    throw std::runtime_error( std::string( "TYPE is missing for the generated column: " )
                                                  .append( std::string( rjCol[ COLUMN_NAME ] ) )
                                                  .append( " in the table " )
                                                  .append( getTableName( ) ) );
                if( rjCol.count( AS ) == 0 )
                    throw std::runtime_error( std::string( "AS is missing for the generated column: " )
                                                  .append( std::string( rjCol[ COLUMN_NAME ] ) )
                                                  .append( " in the table " )
                                                  .append( getTableName( ) ) );
                // append to the SQL statement
                sStmt
                    // Insert separating comma. since the last regular column does not have a comma we can
                    // always safely insert a separating comma
                    .append( ", " )
                    .append( std::string( rjCol[ COLUMN_NAME ] ) ) // column name
                    .append( " " )
                    .append( std::string( rjCol[ TYPE ] ) ) // column type

                    // Generated columns are columns that are computed from other columns
                    // InnoDb supports indices on VIRTUAL generated columns, so there is no need to make use of
                    // the STORED variant. -> actually it doesn't...
                    // The user is supposed to supply an expression that computes the generated column
                    // see:
                    // https://dev.mysql.com/doc/refman/5.7/en/create-table-generated-columns.html
                    // https://dev.mysql.com/doc/refman/5.7/en/create-table-secondary-indexes.html
                    .append( " AS ( " )
                    .append( std::string( rjCol[ AS ] ) )
                    .append( " ) " )
                    .append( rjCol.count( STORED ) == 0 || !rjCol[ STORED ] ? "VIRTUAL" : "STORED" );
                // generated columns can have constraints as well
                if( rjCol.count( CONSTRAINTS ) )
                    // CONSTRAINTS must be plain text describing all constraints
                    sStmt.append( " " ).append( std::string( rjCol[ CONSTRAINTS ] ) );
            } // for

        // Primary Key (for composite form of primary key)
        if( jTableDef.count( PRIMARY_KEY ) )
            sStmt.append( ", PRIMARY KEY (" ).append( std::string( jTableDef[ PRIMARY_KEY ] ) ).append( ")" );

        // Foreign Key definitions
        for( auto& rjItem : jTableDef.items( ) )
            if( rjItem.key( ) == FOREIGN_KEY )
            {
                // Found foreign key definition
                auto& rjVal = rjItem.value( );
                if( rjVal.count( COLUMN_NAME ) != 1 || rjVal.count( REFERENCES ) != 1 )
                    throw std::runtime_error( "JSON table definition for table: " + getTableName( ) +
                                              "\nFOREIGN_KEY definition does not comprise exactly "
                                              "one COLUMN_NAME item as well as REFERENCE item." );
                sStmt.append( ", FOREIGN KEY (" )
                    .append( std::string( rjVal[ COLUMN_NAME ] ) )
                    .append( ")" )
                    .append( " REFERENCES " )
                    .append( std::string( rjVal[ REFERENCES ] ) );
            } // if
        // Close statement
        sStmt.append( ")" );
#ifdef SQL_VERBOSE
        std::cout << "SQL table creation statement: " << sStmt << std::endl;
#endif
        return sStmt;
    } // method


#ifdef POSTGRESQL
    /** @brief Implementation part of getValuesStmt for PG (meta) */
    template <typename... ArgTypes, std::size_t... Is>
    void getValueStmtImpl_PG_( std::string& sStmtText, size_t& uiArgCount, std::index_sequence<Is...> ) const
    {
        ( ( sStmtText.append( ( Is == 0 ? "" : "," ) )
                .append( DBCon::DBImplForward::TypeTranslator::template getPlaceholderForType<ArgTypes>(
                    "$" + std::to_string( uiArgCount++ ) ) ) ),
          ... );
    } // meta

    /** @brief Implementation part of getValuesStmt for PG (core) */
    template <typename FstCol, typename... RemainCols>
    void getValueStmtImpl_PG( std::string& sStmtText, size_t& uiArgCount ) const
    {
        if( std::is_same<FstCol, PGBigSerial>::value )
        {
            // The first column is of type PGBigSerial (indicating a primary key column)
            // PG requires an injection of the keyword 'DEFAULT' for such columns.
            sStmtText.append( "DEFAULT," );
            getValueStmtImpl_PG_<RemainCols...>( sStmtText, uiArgCount, std::index_sequence_for<RemainCols...>{} );
        }
        else
            getValueStmtImpl_PG_<ColTypes...>( sStmtText, uiArgCount, std::index_sequence_for<ColTypes...>{} );
    } // method
#else // MySQL
    /** @brief Implementation part of getValuesStmt */
    template <std::size_t... Is> void getValueStmtImpl_MYSQL( std::string& sStmtText, std::index_sequence<Is...> ) const
    {
        ( ( sStmtText.append( ( Is == 0 ? "" : "," ) )
                .append( DBCon::DBImplForward::TypeTranslator::template getPlaceholderForType<ColTypes>( ) ) ),
          ... );
    } // meta
#endif
    /** @brief Delivers the "VALUES"-part of an insertion statement */
    inline std::string getValuesStmt( size_t& uiArgCount ) const
    {
        std::string sStmtText( "(" );
#ifdef POSTGRESQL
        // Values for columns of type PGBigSerial are automatically inserted by PG
        getValueStmtImpl_PG<ColTypes...>( sStmtText, uiArgCount );
#else // MySQL
        getValueStmtImpl_MYSQL( sStmtText, std::index_sequence_for<ColTypes...>{} );
#endif
        // For INSERT, REPLACE, and UPDATE, if a generated column is inserted into, replaced, or updated
        // explicitly, the only permitted value is DEFAULT. this assumes that the create table statement always
        // appends all generated columns at the end
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
    virtual std::string makeInsertStmt( const size_t uiNumVals = 1 ) const
    {
        // Number of columns
        assert( this->rjTableCols.size( ) == sizeof...( ColTypes ) );

        // DEL const std::string sValuesPartOfStmt( getValuesStmt( ) );
        std::string sStmt = "INSERT INTO ";
        sStmt.append( getTableName( ) ).append( " VALUES " );

        size_t uiArgCount = 1; // used in PostgreSQL for arg counting
        for( size_t uiItrRow = 0; uiItrRow < uiNumVals; )
        {
            sStmt.append( getValuesStmt( uiArgCount ) );

            // insert separating comma
            if( ++uiItrRow < uiNumVals )
                sStmt.append( ", " );
        } // outer for
#if SQL_VERBOSE
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

  public:
    /** @brief Drop the table in DB. (Should only be done by the destructor.) */
    void drop( )
    {
        pDB->execSQL( makeTableDropStmt( ) );
    } // method

  protected:
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

    /** @brief Inserts the argument pack as fresh row into the table without guaranteeing type correctness at
     * compile time.
     *  @detail This form of the insert does not guarantees that there are no problems with column types at
     * runtime. It allows the injection of NULL values into an insertion by passing (void *)NULL at th column
     * position, where the NULL value shall occurs. Don't use insertNonSafe, if you do not need explicit NULL
     * values in the database.
     */
    template <typename... ArgTypes> inline void insertNonSafe( const ArgTypes&... args )
    {
        // static_assert( sizeof...( ArgTypes ) == sizeof...( ColTypes ) );
#ifdef POSTGRESQL
        pInsertStmt->bindAndExec( args... );
#else // only for MySQL
        auto iAffectedRows = pInsertStmt->bindAndExec( args... );
        if( iAffectedRows != 1 )
            throw std::runtime_error( "SQL DB insert. iAffectedRows != 1. (" + std::to_string( iAffectedRows ) + ")" );
#endif
    } // method

    /** @brief: Type-safely insert the argument pack as fresh row into the table. */
    inline SQLTable& insert( const ColTypes&... args )
    {
        this->insertNonSafe( args... );
        return *this;
    } // method

    // Alias for bulk inserter type
    template <size_t BUF_SIZE> using SQLBulkInserterType = SQLBulkInserter<BUF_SIZE, ColTypes...>;

    /** @brief Returns a shared pointer to a bulk inserter for inserting table rows type-safely.
     *   TIP: Definition of a bulk inserter (don't forget the keyword template with GCC!):
     *        auto xBulkInserter = table.template getBulkInserter< size >();
     */
    template <size_t BUF_SIZE> std::shared_ptr<SQLBulkInserterType<BUF_SIZE>> getBulkInserter( )
    {
        return std::make_shared<SQLBulkInserterType<BUF_SIZE>>( *this );
    } // method
#if 0 // File bulk inserters are not used currently
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
#endif
    /** @brief: Delete all rows in the table. */
    SQLTable& deleteAllRows( )
    {
        pDB->execSQL( makeTableDeleteAllRowsStmt( ) );
        return *this;
    } // method

    /** @brief: Creates requested index if not already existing. */
    void addIndex( const json& rjIndexDef )
    {
        SQLIndexView( this, rjIndexDef );
    } // method

    /** @brief: drops requested index if existing. */
    void dropIndex( const json& rjIndexDef )
    {
        SQLIndexDropView( this, rjIndexDef );
    } // method

    /** @brief holds the columns types.
     *  @details
     *  since a parameter pack cannot be held in a using, we hold the type TypePack<ColTypes...>.
     *  pack is a completely empty struct.
     *  @see get_inserter_container_module.h, there this is used to extract the types of all column from a tabletype
     * that is given as a template parameter.
     */
    using ColTypesForw = TypePack<ColTypes...>;

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
    virtual ~SQLTable( )
    {
        do_exception_safe( [&]( ) {
            if( pDB && bDropOnDestr )
                this->drop( );
        } );
    } // virtual destructor (required by clang)

}; // class (SQLTable)

// Used for indicating that the bulk inserter shall reserve a batch of keys.
const bool RESERVE_BATCH_OF_PRIKEY = true;
const bool NO_BATCH_OF_PRIKEY = false;
#if 0
/** @brief Variation of the normal SQL table having an automatic primary key.
 *  @details PriKeyType should be either std::nullptr_t or PriKeyDefaultType. If PriKeyType is
 * PriKeyDefaultType, the primary key is incremented by the SQL library itself. Primary key is always the first
 * column called 'id'.
 */
template <typename DBCon, typename PriKeyType, typename... ColTypes>
class SQLTableWithPriKey : public SQLTable<DBCon, PriKeyType, ColTypes...>
// ORIG: class SQLTableWithPriKey : public SQLTable<DBCon, PriKeyDefaultType, ColTypes...>
{
  private:
    /** @brief Returns if AUTO_INCREMENT is required as part of the primary key definition or not. */
    bool autoIncrRequired( Identity<std::nullptr_t> _ )
    {
        return true;
    } // method

    bool autoIncrRequired( Identity<PriKeyDefaultType> _ )
    {
        return false;
    } // method

    /** @brief Inject primary key column definition into table definition. */
    json inject( const json& rjTableDef )
    {
        // Create deep copy of the table definition
        json jTableDef( rjTableDef );

        // Primary key on first row:
        if( jTableDef.count( TABLE_COLUMNS ) )
            jTableDef[ TABLE_COLUMNS ].insert(
                jTableDef[ TABLE_COLUMNS ].begin( ),
                json::object(
                    {{COLUMN_NAME, "id"},
                     {CONSTRAINTS, std::string( "NOT NULL " ) +
                                       ( autoIncrRequired( Identity<PriKeyType>( ) ) ? " AUTO_INCREMENT " : "" ) +
                                       " UNIQUE PRIMARY KEY"}} ) );
        return jTableDef;
    } // method

    // private attributes
    //-- typename SQLTable<DBCon, PriKeyType, ColTypes...>::SQLIndexView xPriKeyIndexViewType;
    std::shared_ptr<std::atomic<PriKeyDefaultType>>
        pPriKeyCntr; // atomic incrementing counter. (Shared with other instances of the table.)

    /** @brief Insert the argument pack as fresh row into the table.
     *  Only used if the primary key is computed by the library.
     */
    template <typename... ArgTypes>
    inline PriKeyDefaultType dispatchedInsert( Identity<PriKeyDefaultType> _, ArgTypes&&... args )
    {
        // Fetch the current primary key counter and increments it by one
        // fetch_add(1) is a single atomic read-modify-write operation.
        auto uiPriKeyValue = pPriKeyCntr->fetch_add( 1 );
        SQLTable<DBCon, PriKeyDefaultType, ColTypes...>::insertNonSafe( uiPriKeyValue,
                                                                        std::forward<ArgTypes>( args )... );
        return uiPriKeyValue;
    } // method

    /** @brief Insert the argument pack as fresh row into the table, if the primary key is AUTO_INCREMENT. */
    template <typename... ArgTypes>
    inline PriKeyDefaultType dispatchedInsert( Identity<std::nullptr_t> _, ArgTypes&&... args )
    {
        std::lock_guard<std::mutex> xGuard( this->pDBCon->pGlobalInsertLock );
        // Thread-safe region (till end of method).
        SQLTable<DBCon, PriKeyDefaultType, ColTypes...>::insertNonSafe( nullptr, std::forward<ArgTypes>( args )... );
        return this->pDB->getLastAutoIncrVal( );
    } // method

    // private attributes
    std::shared_ptr<DBCon> pDBCon;

  public:
    /** @brief holds the columns types.
     *  @details
     *  since a parameter pack cannot be held in a using, we hold the type TypePack<ColTypes...>.
     *  pack is a completely empty struct.
     *  @see get_inserter_container_module.h, there this is used to extract the types of all column from a
     * tabletype that is given as a template parameter.
     */
    using ColTypesForw = TypePack<ColTypes...>;

    SQLTableWithPriKey( const SQLTableWithPriKey& ) = delete; // no table copies

    /* Constructor */
    SQLTableWithPriKey( std::shared_ptr<DBCon> pDBCon, const json& rjTableDef )
        : SQLTable<DBCon, PriKeyDefaultType, ColTypes...>( pDBCon, inject( rjTableDef ) ),
          pDBCon( pDBCon ) // call superclass constructor
                           // xPriKeyIndexViewType( this, json{ { "INDEX_NAME", "pri_key_idx" }, {
                           // "INDEX_COLUMNS", "id" } } )
    {
        // Request the current maximum of the primary key counter.
        // If the table is empty, the initial primary key counter value is 0.
        // See:
        // https://stackoverflow.com/questions/15475059/how-to-treat-max-of-an-empty-table-as-0-instead-of-null
        SQLQuery<DBCon, PriKeyDefaultType> xQuery( pDBCon,
                                                   "SELECT COALESCE(MAX(id), 0) FROM " + this->getTableName( ) );
        xQuery.execAndFetch( ); // execute the query and fetch the value
        pPriKeyCntr = xSQLDBGlobalSync.getPtrToPriKeyCnt( this->getTableName( ), xQuery.getVal( ) );
    } // constructor

    /** @brief Insert the argument pack as fresh row into the table.
     *	Returns the primary key of the inserted row as return value.
     *  @details A NULL is inserted on positions, where std::nullptr_t occurs in ArgTypes.
     *  NOTICE: This variant is not typesafe, i.e. the correctness of the arguments types can not be guaranteed
     *at compile time.
     */
    template <typename... ArgTypes> inline PriKeyDefaultType insertNonSafe( ArgTypes&&... args )
    {
        return dispatchedInsert( Identity<PriKeyType>( ), std::forward<ArgTypes>( args )... );
    } // method

    /** @brief Type-safely inserts the argument pack as fresh row into the table.
     *	Returns the primary key of the inserted row as return value.
     */
    inline PriKeyDefaultType insert( const ColTypes&... args )
    {
        return dispatchedInsert( Identity<PriKeyType>( ), args... );
    } // method

    /** @brief Alias for bulk inserter type for tables with primary key */
    template <size_t BUF_SIZE>
    using SQLBulkInserterType =
        typename SQLTable<DBCon, PriKeyDefaultType, ColTypes...>::template SQLBulkInserter<BUF_SIZE, PriKeyType,
                                                                                           ColTypes...>;

    /** @brief Delivers a shared pointer to a bulk inserter for inserting table rows type-safely.
     *   TIP: Definition of a bulk inserter (don't forget the keyword 'template' when using GCC!):
     *        auto xBulkInserter = table.template getBulkInserter< size >();
     */
    template <size_t BUF_SIZE>
    std::shared_ptr<SQLBulkInserterType<BUF_SIZE>> getBulkInserter( bool bReserveBatchOfPriKeys = true )
    {
        return std::make_shared<SQLBulkInserterType<BUF_SIZE>>( *this, pPriKeyCntr, bReserveBatchOfPriKeys );
    } // method


    template <typename SQLDBConPoolType, size_t BUF_SIZE>
    using AsyncSQLBulkInserter =
        typename SQLTable<DBCon, PriKeyDefaultType, ColTypes...>::template AsyncSQLBulkInserter<
            SQLDBConPoolType, BUF_SIZE, PriKeyType, ColTypes...>;


    template <size_t BUF_SIZE, typename SQLDBConPoolPtrType>
    std::shared_ptr<AsyncSQLBulkInserter<SQLDBConPoolPtrType, BUF_SIZE>>
    getAsyncBulkInserter( std::shared_ptr<SQLDBConPoolPtrType> pDBConPool, bool bReserveBatchOfPriKeys = true )
    {
        return std::make_shared<AsyncSQLBulkInserter<SQLDBConPoolPtrType, BUF_SIZE>>( pDBConPool, *this, pPriKeyCntr,
                                                                                      bReserveBatchOfPriKeys );
    } // method
}; // class  SQLTableWithPriKey
#endif
/** @brief Variation of the normal SQL table having an automatic primary key.
 *  @details PriKeyType should be either std::nullptr_t or PriKeyDefaultType. If PriKeyType is
 * PriKeyDefaultType, the primary key is incremented by the SQL library itself. Primary key is always the first
 * column called 'id'.
 */
template <typename DBCon, typename PriKeyType, typename... ColTypes>
class SQLTableWithPriKey : public SQLTable<DBCon, PriKeyType, ColTypes...>
{
  private:
    /** @brief Inject primary key column definition into table definition. */
    json inject( const json& rjTableDef, bool autoIncrRequired )
    {
        // Create deep copy of the table definition
        json jTableDef( rjTableDef );

        // Primary key on first row:
        if( jTableDef.count( TABLE_COLUMNS ) )
            jTableDef[ TABLE_COLUMNS ].insert(
                jTableDef[ TABLE_COLUMNS ].begin( ),
                json::object(
                    {{COLUMN_NAME, "id"},
                     {CONSTRAINTS, std::string( "NOT NULL " ) + ( autoIncrRequired ? " AUTO_INCREMENT " : "" ) +
                                       " UNIQUE PRIMARY KEY"}} ) );
        return jTableDef;
    } // method

  public:
    /** @brief holds the columns types.
     *  @details
     *  since a parameter pack cannot be held in a using, we hold the type TypePack<ColTypes...>.
     *  pack is a completely empty struct.
     *  @see get_inserter_container_module.h, there this is used to extract the types of all column from a
     * tabletype that is given as a template parameter.
     */
    using ColTypesForw = TypePack<ColTypes...>;

    SQLTableWithPriKey( const SQLTableWithPriKey& ) = delete; // no table copies

    /* Constructor */
    SQLTableWithPriKey( std::shared_ptr<DBCon> pDBCon, const json& rjTableDef, bool autoIncrRequired )
        : SQLTable<DBCon, PriKeyType, ColTypes...>( pDBCon, inject( rjTableDef, autoIncrRequired ) )
    // call superclass constructor
    {} // constructor
}; // class  SQLTableWithPriKey

template <typename DBCon, typename... ColTypes>
class SQLTableWithLibIncrPriKey : public SQLTableWithPriKey<DBCon, PriKeyDefaultType, ColTypes...>
{
  private:
    // atomic global incrementing counter. (Shared with other instances of the table.)
    std::shared_ptr<std::atomic<PriKeyDefaultType>> pPriKeyCntr;

    /** @brief Insert the argument pack as fresh row into the table.
     *  Only used if the primary key is computed by the library.
     */
    template <typename... ArgTypes> inline PriKeyDefaultType dispatchedInsert( ArgTypes&&... args )
    {
        // Fetch the current primary key counter and increments it by one
        // fetch_add(1) is a single atomic read-modify-write operation.
        auto uiPriKeyValue = pPriKeyCntr->fetch_add( 1 );
        SQLTable<DBCon, PriKeyDefaultType, ColTypes...>::insertNonSafe( uiPriKeyValue,
                                                                        std::forward<ArgTypes>( args )... );
        return uiPriKeyValue;
    } // method

  public:
    /* Constructor */
    SQLTableWithLibIncrPriKey( std::shared_ptr<DBCon> pDBCon, const json& rjTableDef )
        : SQLTableWithPriKey<DBCon, PriKeyDefaultType, ColTypes...>( pDBCon, rjTableDef, false /* NO AUTO_INCR */ )
    // call superclass constructor
    {
        // Request the current maximum of the primary key counter.
        // If the table is empty, the initial primary key counter value is 0.
        // See:
        // https://stackoverflow.com/questions/15475059/how-to-treat-max-of-an-empty-table-as-0-instead-of-null
        SQLQuery<DBCon, PriKeyDefaultType> xQuery( pDBCon,
                                                   "SELECT COALESCE(MAX(id), 0) FROM " + this->getTableName( ) );
        // xQuery.execAndFetch( ); // execute the query and fetch the value
        auto uiCurrPriKey = xQuery.execAndGetValue( );
        pPriKeyCntr = xSQLDBGlobalSync.getPtrToPriKeyCnt( this->getTableName( ), uiCurrPriKey /* xQuery.getVal( ) */ );
    } // constructor

    /** @brief Insert the argument pack as fresh row into the table with NULL values being possible.
     *	Returns the primary key of the inserted row as return value.
     *  @details A NULL is inserted on positions, where std::nullptr_t occurs in ArgTypes.
     *  NOTICE: This variant is not typesafe, i.e. the correctness of the arguments types can not be guaranteed
     *  at compile time.
     */
    template <typename... ArgTypes> inline PriKeyDefaultType insertNonSafe( ArgTypes&&... args )
    {
        return dispatchedInsert( std::forward<ArgTypes>( args )... );
    } // method

    /** @brief Type-safely inserts the argument pack as fresh row into the table.
     *	Returns the primary key of the inserted row as return value.
     */
    inline PriKeyDefaultType insert( const ColTypes&... args )
    {
        return dispatchedInsert( args... );
    } // method

    /** @brief Alias for bulk inserter type for tables with primary key */
    template <size_t BUF_SIZE>
    using SQLBulkInserterType =
        typename SQLTable<DBCon, PriKeyDefaultType, ColTypes...>::template SQLBulkInserter<BUF_SIZE, PriKeyDefaultType,
                                                                                           ColTypes...>;

    /** @brief Delivers a shared pointer to a bulk inserter for inserting table rows type-safely.
     *   TIP: Definition of a bulk inserter (don't forget the keyword 'template' when using GCC!):
     *        auto xBulkInserter = table.template getBulkInserter< size >();
     */
    template <size_t BUF_SIZE>
    std::shared_ptr<SQLBulkInserterType<BUF_SIZE>> getBulkInserter( bool bReserveBatchOfPriKeys = true )
    {
        return std::make_shared<SQLBulkInserterType<BUF_SIZE>>( *this, pPriKeyCntr, bReserveBatchOfPriKeys );
    } // method

    template <typename SQLDBConPoolType, size_t BUF_SIZE>
    using AsyncSQLBulkInserter =
        typename SQLTable<DBCon, PriKeyDefaultType, ColTypes...>::template AsyncSQLBulkInserter<
            SQLDBConPoolType, BUF_SIZE, PriKeyDefaultType, ColTypes...>;

    template <size_t BUF_SIZE, typename SQLDBConPoolPtrType>
    std::shared_ptr<AsyncSQLBulkInserter<SQLDBConPoolPtrType, BUF_SIZE>>
    getAsyncBulkInserter( std::shared_ptr<SQLDBConPoolPtrType> pDBConPool, bool bReserveBatchOfPriKeys = true )
    {
        return std::make_shared<AsyncSQLBulkInserter<SQLDBConPoolPtrType, BUF_SIZE>>( pDBConPool, *this, pPriKeyCntr,
                                                                                      bReserveBatchOfPriKeys );
    } // method
}; // class

#ifdef POSTGRESQL
using AutoIncrPriKeyTypeForTblDef = PGBigSerial;
using AutoIncrPriKeyTypeForBulkIns = PGBigSerial;
#else
using AutoIncrPriKeyTypeForTblDef = PriKeyDefaultType;
using AutoIncrPriKeyTypeForBulkIns = std::nullptr_t; /* With AUTO_INCR */
#endif

template <typename DBCon, typename... ColTypes>
class SQLTableWithAutoPriKey : public SQLTableWithPriKey<DBCon, AutoIncrPriKeyTypeForTblDef, ColTypes...>
{
  private:
    /** @brief Insert the argument pack as fresh row into the table, if the primary key is AUTO_INCREMENT.
     * Returns the value of the primary key of the freshly inserted row.
     */
    template <typename... ArgTypes> inline PriKeyDefaultType insert__( ArgTypes&&... args )
    {
        // This lock might not be necessary with POSTGRESQL
        std::lock_guard<std::mutex> xGuard( this->pDBCon->pGlobalInsertLock );
        // Thread-safe region (till end of method).
        SQLTable<DBCon, AutoIncrPriKeyTypeForTblDef, ColTypes...>::insertNonSafe(
#ifndef POSTGRESQL // MYSQL
            nullptr, // In MySQL we must additionally pass a nullptr here for the primary key.
#endif
            std::forward<ArgTypes>( args )... );
#ifdef POSTGRESQL
        return xCurrvalQuery.execAndGetValue( );
#else
        return this->pDB->getLastAutoIncrVal( );
#endif
    } // method

    std::shared_ptr<DBCon> pDBCon;

#ifdef POSTGRESQL
    // See: https://stackoverflow.com/questions/2944297/postgresql-function-for-last-inserted-id
    // Query for requesting the
    SQLQuery<DBCon, int64_t> xCurrvalQuery;
#endif

  public:
    /* Constructor */
    SQLTableWithAutoPriKey( std::shared_ptr<DBCon> pDBCon, const json& rjTableDef )
        : SQLTableWithPriKey<DBCon, AutoIncrPriKeyTypeForTblDef, ColTypes...>( pDBCon, rjTableDef,
#ifdef POSTGRESQL
                                                                               false /* NO AUTO_INCR */
#else
                                                                               true /* With AUTO_INCR */
#endif
                                                                               ),
          pDBCon( pDBCon )
#ifdef POSTGRESQL
          ,
          xCurrvalQuery( pDBCon, "SELECT currval(pg_get_serial_sequence('" + this->getTableName( ) + "','id'))" )
#endif
    // call superclass constructor
    {} // constructor

    /** @brief Insert the argument pack as fresh row into the table.
     *	Returns the primary key of the inserted row as return value.
     *  @details A NULL is inserted on positions, where std::nullptr_t occurs in ArgTypes.
     *  NOTICE: This variant is not typesafe, i.e. the correctness of the arguments types can not be guaranteed
     *at compile time.
     */
    template <typename... ArgTypes> inline PriKeyDefaultType insertNonSafe( ArgTypes&&... args )
    {
        return insert__( std::forward<ArgTypes>( args )... );
    } // method

    /** @brief Type-safely inserts the argument pack as fresh row into the table.
     *	Returns the primary key of the inserted row as return value.
     */
    inline PriKeyDefaultType insert( const ColTypes&... args )
    {
#ifdef POSTGRESQL
        return insert__( args... );
#else // MYSQL
        return insert__( /* nullptr, */ args... ); // MySQL requires the explicit insertion of the nullptr
#endif
    } // method

    /** @brief Alias for bulk inserter type for tables with primary key */
    template <size_t BUF_SIZE>
    using SQLBulkInserterType =
        typename SQLTable<DBCon, AutoIncrPriKeyTypeForTblDef,
                          ColTypes...>::template SQLBulkInserter<BUF_SIZE, AutoIncrPriKeyTypeForBulkIns, ColTypes...>;

    /** @brief Delivers a shared pointer to a bulk inserter for inserting table rows type-safely.
     *   TIP: Definition of a bulk inserter (don't forget the keyword 'template' when using GCC!):
     *        auto xBulkInserter = table.template getBulkInserter< size >();
     */
    template <size_t BUF_SIZE>
    std::shared_ptr<SQLBulkInserterType<BUF_SIZE>> getBulkInserter( bool bReserveBatchOfPriKeys = true )
    {
        // #ifdef POSTGRESQL
        //         return std::make_shared<SQLBulkInserterType<BUF_SIZE>>( *this, nullptr, bReserveBatchOfPriKeys );
        // #else // MySQL
        return std::make_shared<SQLBulkInserterType<BUF_SIZE>>( *this, nullptr, bReserveBatchOfPriKeys );
        // #endif
    } // method
}; // class

/** @brief Table with primary key of type PriKeyDefaultType, where the incrementing is done by the library.
 *  This kind of tables allows the bulk-inserter to return the primary key immediately after insertion.
 */
// template <typename DBCon, typename... ColTypes>
// using SQLTableWithLibIncrPriKey = SQLTableWithPriKey<DBCon, PriKeyDefaultType, ColTypes...>;

// template <typename DBCon, typename... ColTypes>
// using SQLTableWithAutoPriKey = SQLTableWithPriKey<DBCon, std::nullptr_t, ColTypes...>;

/** @brief Informs about the database in general */
template <typename DBConType> class SQLDBInformer
{
  private:
    // private attributes
    std::shared_ptr<DBConType> pDBCon; // DB connection used by SQLDBInformer

  public:
    SQLDBInformer( std::shared_ptr<DBConType> pDBCon ) : pDBCon( pDBCon )
    {} // constructor

    /** @brief Returns a vector of strings holding all schema-names */
    std::vector<std::string> getAllSchemas( )
    {
        SQLQuery<DBConType, std::string> xQuery( pDBCon, "SHOW SCHEMAS" );
        return xQuery.template executeAndStoreInVector<0>( );
    } // method
}; // class SQLDBInformer

// template <typename DBConType> using SQLDBInformer = _SQLDBInformer<std::shared_ptr<DBConType>>;

#define DEFAULT_DBNAME "defaultDB"

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

        /** @brief Constructs a transaction. */
        GuardedTransaction( SQLDB<DBImpl>& rHost, std::shared_ptr<bool> pHostTombStone )
            : rHost( rHost ), pHostTombStone( pHostTombStone ), bStateCommitted( true )
        {
            this->start( );
        } // constructor

        /** @brief Commit the transaction */
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

        /** @brief Start the transaction */
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

        /** @brief Commits and destructs transaction. */
        ~GuardedTransaction( )
        {
            this->commit( );
        } // destructor
    }; // class (GuardedTransaction)

    std::shared_ptr<bool> pTombStone; // So long the database connection is alive, the stone says false.
                                      // As soon as the DB is gone, the stone says true.

    /** @brief Extracts the schema name from JSON definition of the connection. Uses DEFAULT_DBNAME as schema
     * name, if the JSON does not define any name.
     */
    std::string getSchemaName( const json& jDBConfig )
    {
        if( jDBConfig.count( SCHEMA ) && jDBConfig[ SCHEMA ].count( NAME ) )
            return jDBConfig[ SCHEMA ][ NAME ].get<std::string>( );
        return DEFAULT_DBNAME;
    } // method

    /** @brief Extracts the value of the DROP_ON_CLOSURE flag, if existing. Otherwise returns false. */
    const std::string _FLAGS = "FLAGS"; // @todo fix flags
    bool getDropOnClosureFlag( const json& jDBConfig )
    {
        if( jDBConfig.count( SCHEMA ) && jDBConfig[ SCHEMA ].count( _FLAGS ) )
            for( auto rsString : jDBConfig[ SCHEMA ][ _FLAGS ] )
                if( rsString == DROP_ON_CLOSURE )
                    return true;
        return false;
    } // method

  public:
    typedef std::unique_ptr<GuardedTransaction> uniqueGuardedTrxnType;
    typedef std::shared_ptr<GuardedTransaction> sharedGuardedTrxnType;
    typedef DBImpl DBImplForward; // Forwarding of the template parameter type
    const std::string sConId; // unique connection ID with respect to the current machine (used by the FileBulkInserter)
    const std::string sSchemaName; // the schema name used by the connection
    const bool bDropOnClosure; // this is true the schema shall self delete in its destructor

    SQLDB( const SQLDB& ) = delete; // no DB connection copies

    // For serializing all table insertions with respect to the current DB-connection.
    // Necessary for getting the correct last insertion ID in case of an auto-increment column. (e.g. with
    // MySQL)
    // FIXME: Remove the lock, because it is not necessary.
    std::mutex pGlobalInsertLock;

    SQLDB( const json& jDBConData = json{} )
        : DBImpl( jDBConData ),
          pTombStone( std::make_shared<bool>( false ) ), // initialize tombstone
          sConId( intToHex( reinterpret_cast<uint64_t>( this ) ) ), // use the address for id creation
          sSchemaName( getSchemaName( jDBConData ) ),
          bDropOnClosure( getDropOnClosureFlag( jDBConData ) )
    {
        // std::cout << jDBConData << std::endl;
        // Register the selected schema with global warden and select it for use.
        xSQLDBGlobalSync.registerSchema( sSchemaName );
        DBImpl::useSchema( sSchemaName );
    } // constructor

    /** @brief Initialize DB connection using a given schema name.
     *  @details This constructor is exported to python.
     */
    SQLDB( std::string sSchemaName )
        : SQLDB( json{{SCHEMA, {{NAME, sSchemaName}}},
#ifdef USE_PG
                      { CONNECTION,
                        { {HOSTNAME, "localhost"},
                          {USER, "postgres"},
                          {PASSWORD, "admin"},
                          { PORT,
                            0 } } }
#endif
#ifdef USE_MSQL
                      { CONNECTION,
                        { {HOSTNAME, "localhost"},
                          {USER, "root"},
                          {PASSWORD, "admin"},
                          { PORT,
                            0 } } }
#endif
          } )
    {}

    /** @brief Destructs connection in an exception safe way ... */
    ~SQLDB( )
    {
        // Throwing an exception in a destructor results in undefined behavior.
        doNoExcept( [this] {
            // Unregister the schema in the global visor
            auto uiRemainCons = xSQLDBGlobalSync.unregisterSchema( sSchemaName );

            // If I am the last connection using my schema and the schema shall be dropped finally,
            // it is my job to do it now ...
            if( ( uiRemainCons == 0 ) && this->bDropOnClosure )
            {
                xSQLDBGlobalSync.doSynchronized( [this] {
                    // DEBUG: std::cout << "Do schema drop synchronized in Destructor" << std::endl;
                    this->dropSchema( this->sSchemaName );
                } ); // doSynchronized
            } // if
        } );
        // Inform about the "death" of the database connection.
        *pTombStone = true;
        // Now the destructor of the DB implementor is called ...
    }; // destructor

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

    /** @brief Drop Database Schema.
     *  @details
     *  @note Used id with care!
     */
    void dropSchema( const std::string& rsSchemaName )
    {
#ifdef POSTGRESQL
        DBImpl::execSQL( "DROP SCHEMA " + rsSchemaName + " CASCADE" );
#else
        DBImpl::execSQL( "DROP DATABASE " + rsSchemaName );
#endif
    } // method

    /** @brief This function should be redefined in pooled SQL connections so that the function is executed in
     * mutex protected environment. In a single connection we simply execute the function.
     *  @todo Get rid of doPoolSafe here by doing the locking via the master sync object.
     */
    template <typename F> void doPoolSafe( F&& func )
    {
        // DEBUG: std::cout << "doPoolSafe in SQLDB ..." << std::endl;
        func( );
    } // method
}; // SQLDB
