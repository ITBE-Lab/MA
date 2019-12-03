/* This file is part of MA.
 * Authors: Markus Schmidt and Arne Kutzner
 * Created: Oct. 2018
 * MIT License
 * @file sqlite3.h
 */
#define USE_NEW_DB_CODE
#ifdef USE_NEW_DB_CODE
#include <db_sql.h>
#else

#pragma once

#include "debug.h"
#include <CppSQLite3.h>
#include <cassert>
#include <iostream>
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <vector>

/* Forward declaration of class for SQL statements
 */
class CppSQLiteExtStatementParent;

template <class... Types> class CppSQLiteExtQueryStatementParent;

template <class... Types> struct SQLQueryTemplateString;

/* Enumerated type for transporting data with respect to the database opening mode.
 */
enum enumSQLite3DBOpenMode
{
    eOPEN_DB,
    eCREATE_DB
}; // enum

/* METAPROGRAMMING
 * Unpack a tuple to call a matching function.
 * FIX ME: Better naming scheme, more documentation.
 * http://stackoverflow.com/questions/7858817/unpacking-a-tuple-to-call-a-matching-function-pointer
 * FIXME: Use the construct from the standard
 */
template <int...> struct Seq
{};

template <int N, int... S> struct Gens : Gens<N - 1, N - 1, S...>
{};

template <int... S> struct Gens<0, S...>
{
    typedef Seq<S...> type;
};

/* METAPROGRAMMING
 * For a solution without std::function look over here: // DEPRECATED - use the below form now.
 * http://stackoverflow.com/questions/9535680/functions-functors-as-template-parameters-can-they-be-stored
 * FIXME: Use the construct from the standard
 */
template <typename... Args> struct TupleUnpackToParameterByReference
{ /* CHECK ME: Shouldn't we use a reference here? Hmm.. There is some discussion about on stackoverflow.
   */
    const std::function<void( const Args&... )>& func;

    template <int... S> void callFunc( const std::tuple<Args...>& params, Seq<S...> )
    {
        func( std::get<S>( params )... );
    }

    void operator( )( const std::tuple<Args...>& params )
    {
        callFunc( params, typename Gens<sizeof...( Args )>::type( ) );
    }

    TupleUnpackToParameterByReference( const std::function<void( const Args&... )>& func ) : func( func )
    {} // constructor
}; // struct

/* METAPROGRAMMING
 * Unpacks a tuple for function application.
 */
template <typename... Args> struct TupleUnpackAndCallFunctor
{
    const std::tuple<Args...>& _params;

    template <typename Functor, int... S> void callFunc( Functor&& f, Seq<S...> )
    {
        f( std::get<S>( _params )... );
    } // method

    template <typename Functor> void operator( )( Functor&& f )
    {
        callFunc( std::forward<Functor>( f ), typename Gens<sizeof...( Args )>::type( ) );
    } // operator ()

    TupleUnpackAndCallFunctor( const std::tuple<Args...>& params ) : _params( params )
    {} // constructor

    /* Constructor that involves directly the functor call. */
    template <typename Functor>
    TupleUnpackAndCallFunctor( Functor&& f, const std::tuple<Args...>& params ) : _params( params )
    {
        callFunc( std::forward<Functor>( f ), typename Gens<sizeof...( Args )>::type( ) );
    } // constructor
}; // struct

/* METAPROGRAMMING
 * 1. Application of some function f to all elements of some tuple.
 * Iteration starts with index 0.
 * Works with the empty tuple as well!
 * Taken from http://pastebin.com/h0Je0453
 */
template <int I, int TYPE_SIZE, typename Tuple>
struct iterate_over_tuple_impl : public iterate_over_tuple_impl<I + 1, TYPE_SIZE, Tuple>
{
    typedef typename std::tuple_element<I, Tuple>::type tp;

    template <typename Function> void operator( )( Function&& f, Tuple& t )
    {
        /* Application of the function to the i-th element of the tuple.
         * Improvement: Deliver the I and TYPE_SZIE as integral template parameter.
         */
        f( std::get<I>( t ), I, TYPE_SIZE );

        iterate_over_tuple_impl<I + 1, TYPE_SIZE, Tuple>::operator( )( std::forward<Function>( f ), t );
    } // operator
}; // struct

/* METAPROGRAMMING */
template <int I, typename Tuple> struct iterate_over_tuple_impl<I, I, Tuple>
{
    template <typename Function> void operator( )( Function&& f, Tuple& t )
    {}
}; // struct

/* METAPROGRAMMING
 * Fills the element of the given tuple t by repeatedly calling the function f.
 * Starts with tuple element 0 and wok up to tuple element sizeof...(TupleTypes).
 */
template <typename Functor, typename... TupleTypes>
void iterateOverTuple( Functor&& functor, std::tuple<TupleTypes...>& tuple )
{
    iterate_over_tuple_impl<0, sizeof...( TupleTypes ), std::tuple<TupleTypes...>>( )( std::forward<Functor>( functor ),
                                                                                       tuple ); // function call
}; // struct

/* METAPROGRAMMING
 * Like iterateOverTuple but with an additional currying of 2 arguments.
 * (So, iterateOverTupleCurry2 requires a functor that expects two arguments.)
 */
template <typename Functor, typename... TupleTypes>
void iterateOverTupleCurry2( Functor&& functor, std::tuple<TupleTypes...>& tuple )
{
    iterate_over_tuple_impl<0, sizeof...( TupleTypes ), std::tuple<TupleTypes...>>( )(
        std::bind( std::forward<Functor>( functor ), std::placeholders::_1, std::placeholders::_2 ),
        tuple ); // function call
}; // struct

struct GetTupleElement
{
    CppSQLite3Query& rxQuery;

    GetTupleElement( CppSQLite3Query& rxQuery ) : rxQuery( rxQuery )
    {}

    /* int CppSQLite3Query::getIntField(int nField, int nNullValue
     * const char* CppSQLite3Query::fieldValue(int nField) [delivers the field value generally as string]
     * long long CppSQLite3Query::getInt64Field(int nField, long long nNullValue=0)
     * double CppSQLite3Query::getFloatField(int nField, double fNullValue=0.0)
     * const char* CppSQLite3Query::getStringField(int nField, const char* szNullValue="")
     * const unsigned char* CppSQLite3Query::getBlobField(int nField, int& nLen)
     */
    void getColumnElement( int& rValue, const int iFieldIndex )
    {
        rValue = rxQuery.getIntField( iFieldIndex );
    } // method

/* The long datatype has different sizes depending on different platforms.
 * So, we have to carefully distinguish over here.
 */
#ifdef __LP64__
    static_assert( sizeof( long ) == 8, "size of long is not 8 bytes for this compiler" ); // GCC 64 Bit
#else
    static_assert( sizeof( long ) == 4, "size of long is not 4 bytes for this compiler" ); // MS VC++ 32 Bit
#endif
    void getColumnElement( long& rValue, const int iFieldIndex )
    {
#ifdef __LP64__
        rValue = rxQuery.getInt64Field( iFieldIndex ); // GCC 64 Bit
#else
        rValue = rxQuery.getIntField( iFieldIndex ); // MS VC++ 32 Bit
#endif
    } // method

    void getColumnElement( long long& rValue, const int iFieldIndex )
    {
        rValue = rxQuery.getInt64Field( iFieldIndex );
    } // method

    void getColumnElement( unsigned int& rValue, const int iFieldIndex )
    {
        rValue = (unsigned int)rxQuery.getInt64Field( iFieldIndex );
    } // method

    void getColumnElement( double& rValue, const int iFieldIndex )
    {
        rValue = rxQuery.getFloatField( iFieldIndex );
    } // method

    void getColumnElement( std::string& rValue, const int iFieldIndex )
    {
        /* This triggers a copy of the string referred by the received reference
         */
        rValue = rxQuery.getStringField( iFieldIndex );
    } // method

    /* Use this overloading with care because it does not create a copy of the referenced string.
     * The string is available, so long we do not switch to the next row. (see SQLite documentation)
     */
    void getColumnElement( const char*& rValue, const int iFieldIndex )
    {
        rValue = rxQuery.getStringField( iFieldIndex );
    } // method

    /* We extract the first character of the string delivered by SQLite. */
    void getColumnElement( char& rValue, const int iFieldIndex )
    {
        rValue = *rxQuery.getStringField( iFieldIndex );
    } // method

    /* We extract the boolean delivered by SQLite. */
    void getColumnElement( bool& bValue, const int iFieldIndex )
    {
        bValue = rxQuery.getIntField( iFieldIndex ) != 0;
    } // method

    void getColumnElement( SQL_BLOB& rValue, const int iFieldIndex )
    {
        /* Copies the blob content into a vector of bytes */
        int iBlobSize;
        const unsigned char* pBlobRef = rxQuery.getBlobField( iFieldIndex, iBlobSize );
        rValue.fromBlob( pBlobRef, iBlobSize );
    } // method

    /* The trick within the operator is the overloading of getColumnElement.
     * The compiler will make already on compile-time the decision regarding the appropriate function call.
     */
    template <typename TP> void operator( )( TP& rElement, const int I )
    {
        getColumnElement( rElement, I );
    } // operator
}; // struct

struct PrintTupleElement
{
    template <typename TP> void operator( )( TP& rElement )
    {
        getColumnElement( rElement );
    } // operator
}; // struct


/* TODO: here we could work with precompiled statements, in order to boost execution time.
 * Idea: We compile the statement with its first execution.
 *       We use parameter binding in order to get the values into the query(request).
 */
struct SQLRequest
{
    /* The request text itself. */
    const char* pcText;

    /* Internal text-buffer used for the storage of compiles statements. */
    CppSQLite3Buffer textBuffer;

    template <class... ArgTypes> const char* insertArguments( ArgTypes... args )
    {
        return textBuffer.format( pcText, args... );
    } // generic method
}; // struct


/* In this approach we give the types of all columns
 */
template <class... Types> class CppSQLite3IteratorTable
{
    /* We give CppSQLiteDBExtended access to the internal table */
    friend class CppSQLiteDBExtended;

  private:
    /* Here we store the internal table
     */
    std::vector<std::tuple<Types...>> xInternalTable;

  public:
    /* Prevent object copies and assignment for this class. */
    CppSQLite3IteratorTable( CppSQLite3IteratorTable const& ) = delete;
    CppSQLite3IteratorTable& operator=( CppSQLite3IteratorTable const& ) = delete;

    /* The type of a single row in the table.
     */
    typedef std::tuple<Types...> tpRowType;

    template <int C> struct getTuple
    {
        typedef typename std::tuple_element<C, std::tuple<Types...>>::type result_type;

        result_type& operator( )( std::tuple<Types...>& t ) const
        {
            return std::get<C>( t );
        } // method
    }; // struct

    /* Moves the content of the argument table into the current table.
     * Uses move semantics, so the argument table stays in some undefined state.
     */
    void moveContentFromArgumentTableToCurrentTable( CppSQLite3IteratorTable& rxTableRef )
    {
        xInternalTable.reserve( xInternalTable.size( ) + rxTableRef.xInternalTable.size( ) );
        std::move( rxTableRef.xInternalTable.begin( ),
                   rxTableRef.xInternalTable.end( ),
                   std::back_inserter( xInternalTable ) );
    }; // method

    /* Copies the content of column C into a set. */
    template <int C> std::set<typename std::tuple_element<C, std::tuple<Types...>>::type> copyColumnIntoSet( )
    {
        typedef typename std::tuple_element<C, std::tuple<Types...>>::type element_type;

        std::set<element_type> xResultSet;

        for( auto xRow : xInternalTable )
        {
            xResultSet.insert( getTuple<C>( xRow ) );
        } // for

        return xResultSet;
    } // method

    /* Moves the content of column C into a vector. The column shouldn't be used after the moveIntoSet operation. */
    template <int C> std::vector<typename std::tuple_element<C, std::tuple<Types...>>::type> moveColumnIntoVector( )
    {
        typedef typename std::tuple_element<C, std::tuple<Types...>>::type element_type;

        std::vector<element_type> xResultSet;

        for( auto xRow : xInternalTable )
        {
            xResultSet.emplace_back( std::move( getTuple<C>( xRow ) ) );
        } // for

        return xResultSet;
    } // method

    /* Delivers the first column in form of some vector.
     */
    // std::vector< std::tuple_element< 0, std::tuple<Types...> >::type >

    auto begin( ) -> decltype( xInternalTable.begin( ) )
    {
        return xInternalTable.begin( );
    } // method

    auto end( ) -> decltype( xInternalTable.begin( ) )
    {
        return xInternalTable.end( );
    } // method

    auto size( ) -> decltype( xInternalTable.size( ) )
    {
        return xInternalTable.size( );
    } // method

    auto at( size_t uiPosition ) -> decltype( xInternalTable.at( uiPosition ) )
    {
        return xInternalTable.at( uiPosition );
    } // method

    auto firstAt( size_t uiPosition ) -> decltype( std::get<0>( xInternalTable.at( uiPosition ) ) )
    {
        return std::get<0>( xInternalTable.at( uiPosition ) );
    } // method

    /* Constructor
     */
    CppSQLite3IteratorTable( )
    {} // constructor

    /* Destructor
     */
    ~CppSQLite3IteratorTable( )
    {
        vFinalize( );
    } // destructor

  private:
    void vFinalize( )
    {} // method
}; // class

/* Extended from of the basic SQLite DB wrapper.
 */
class CppSQLiteDBExtended : public CppSQLite3DB
{
  private:
    /* If we deliver an already existing table reference as second parameter the method works in an append-mode.
     * I.e. all fresh table rows are pushed to the end of some existing table.
     * Be careful that types does not comprise (const char *), because this results in major trouble.
     */
    template <class... Types>
    std::unique_ptr<CppSQLite3IteratorTable<Types...>> ExecuteQueryTableCore( CppSQLite3Query& xQuery )
    {
        /* We check whether the number of columns is ok. */
        if( xQuery.numFields( ) != ( sizeof...( Types ) ) )
        {
            throw CppSQLite3Exception( CPPSQLITE_ERROR,
                                       "actual number of columns does not match number of types of template",
                                       false // => DONT_DELETE_MSG
            );
        } // if

        /* Creation of the result table */
        std::unique_ptr<CppSQLite3IteratorTable<Types...>> pResultTableRef( new CppSQLite3IteratorTable<Types...> );

        /* We create a functor object for the following repeated iterations over the tuple. */
        auto xGetElementFunctor = GetTupleElement( xQuery );

        /* We iterate over all rows of our table. */
        while( !xQuery.eof( ) )
        {
            { /* We store single rows as tuples. By means of tuples we can keep all type info. */
                pResultTableRef->xInternalTable.emplace_back( std::tuple<Types...>( ) );

                std::tuple<Types...>& xSingleRowAsTuple = pResultTableRef->xInternalTable.back( );

                iterateOverTupleCurry2( xGetElementFunctor, xSingleRowAsTuple );
            } // inner block

            xQuery.nextRow( );
        } // while

        return pResultTableRef;
    } // private method

    /* The start of an immediate transaction can boost the execution speed.
     * We have to work with unique pointers over here, because the members of the type CppSQLiteExtStatement are unknown
     * over here.
     */
    std::unique_ptr<CppSQLiteExtStatementParent> xStatementBeginTransaction = nullptr;
    std::unique_ptr<CppSQLiteExtStatementParent> xStatementEndTransaction = nullptr;

    /* Externally defined method that initialized the above attributes.
     */
    void vInititializeTransactionStatement( void );

  public:
    friend class CppSQLiteExtImmediateTransactionContext;

    /* Here we communicate the construction mode to other objects.
     */
    const enumSQLite3DBOpenMode eDatabaseOpeningMode;

    /* Sophisticated constructor that deletes an existing physical database file before creating a new database.
     * The finally open command creates a fresh database.
     */
    CppSQLiteDBExtended(
        const std::string& rsWorkingDirectory, // directory location of the database
        const std::string& sDataBaseName, // name of the database
        const enumSQLite3DBOpenMode eDatabaseOpeningMode // how to open the database (creation or simply opening)
        )
        : CppSQLite3DB( ), // call superclass constructor
          eDatabaseOpeningMode( eDatabaseOpeningMode )
    {
        /* Construct the full database name.*/
        std::string xDatabaseFileName = rsWorkingDirectory + sDataBaseName;

        /* In the case of a fresh creation delete any already physically existing database file.
         * FIXME: re-implement this without boost...
         */
#if 0
        if( ( eDatabaseOpeningMode == eCREATE_DB ) && ( boost::filesystem::exists( xDatabaseFileName ) ) )
        { /* Rename and remove the existing database file.
           */
            BOOST_LOG_TRIVIAL( info ) << "Remove existing SQLite database file " << xDatabaseFileName;

            boost::filesystem::path xDeleteName = boost::filesystem::path( rsWorkingDirectory ) /= "temporary.delete";
            if( std::rename( xDatabaseFileName.string( ).c_str( ), xDeleteName.string( ).c_str( ) ) )
            {
                std::remove( xDeleteName.string( ).c_str( ) );
            } // if
        } // if
#endif

        /* Finally create/open the database.
         * If the database does not exist, it is created over here.
         */
        this->open( xDatabaseFileName.c_str( ) );

        /* Initialize the transaction statements.
         */
        vInititializeTransactionStatement( );
    } // constructor

    /* Simplified constructor that does not require expect the specification of a working directory.
     */
    CppSQLiteDBExtended( const std::string& sDataBaseName, // name of the database
                         const enumSQLite3DBOpenMode eDatabaseOpeningMode // how to open the database
                         )
        : CppSQLiteDBExtended( "", sDataBaseName, eDatabaseOpeningMode )
    {} // constructor

    /* We need some virtual destructor due to the inheritance
     */
    virtual ~CppSQLiteDBExtended( )
    {} // destructor

    /* TODO: Insert data-type check */
    void vCreateTableFromTextFile(
        const std::vector<std::string>& xDatabaseColumns, // The names of the database columns
        const char cDelimiter,
        const char* sTableName,
        const std::function<bool( const char* )>& fCallbackComment, // decides whether we have a comment
        const std::function<bool( const char* )>&
            fCallbackRepresentsNull, // decides whether a column element shall be interpreted as NULL
        const std::vector<std::pair<std::string, std::string>>& rxFileNameVector, // (filename, content of first column)
        bool bInsertColumnIdAsPrimaryKey = true // automatically we insert an id row as primary key
    ); // method prototype

    /* Plain table creation */
    void vCreateTable( const std::vector<std::string>& xDatabaseColumns, // The names of the database columns
                       const char* sTableName, // Database same
                       bool bInsertColumnIdAsPrimaryKey = true,
                       const std::vector<std::string>& vConstraints = {},
                       const bool bWithoutRowId = false ); // method prototype

    /* We query a result table.
     * The parameter types correspond to the types of the columns.
     * The statement is currently still given as text, later we should work with pre-compiled statements here.
     */
    template <class... Types>
    std::unique_ptr<CppSQLite3IteratorTable<Types...>> queryTable( const char* pcSQLStatement )
    {
        /* We create a SQL query object by compiling the statement.
         */
        CppSQLite3Query xQuery = execQuery( pcSQLStatement );

        return ExecuteQueryTableCore<Types...>( xQuery );
    } // method

    /* Here we can infer the types directly from the SQL Request.
     * We work with some fresh copy of the request object. (Quite inefficient)
     */
    template <class... Types, class... ArgTypes>
    std::unique_ptr<CppSQLite3IteratorTable<Types...>> queryTable( SQLRequest request, const ArgTypes&... args )
    {
        /* We forward the query to the code that works with requests in the form of plain strings.
         */
        return queryTable<Types...>( request.insertArguments( args... ) );
    } // public method

    /* Query statement for some pre-compiled query statement, incl. all argument binding.
     * Here we need an r-value reference, because
     *  - const would not be really ok with CppSQLiteExtQueryStatement.
     *  - the construction of universal references is not possible.
     * If pResultTableRef is not a null-pointer the query works in an append-mode
     */
    template <class... Types, class... ArgTypes>
    std::unique_ptr<CppSQLite3IteratorTable<Types...>>
    getTable( CppSQLiteExtQueryStatementParent<Types...>&& rCompiledQuery, ArgTypes&&... args )
    {
        /* We bind all arguments and execute the statement.
         * This process can fail and throw an exception.
         */
        CppSQLite3Query xQuery = rCompiledQuery.bindAndExecQuery( std::forward<ArgTypes>( args )... );
        return ExecuteQueryTableCore<Types...>( xQuery );
    } // public method

    /* Work around. Necessary due to missing const's in CppSQLiteQueryStatement.
     * (pre-complied, but with l-value reference instead of r-value reference)
     */
    template <class... Types, class... ArgTypes>
    std::unique_ptr<CppSQLite3IteratorTable<Types...>>
    getTable( CppSQLiteExtQueryStatementParent<Types...>& rCompiledQuery, ArgTypes&&... args )
    {
        return getTable( std::move( rCompiledQuery ), std::forward<ArgTypes>( args )... );
    } // public method

    /* Query on the foundation of a one time use of some template.
     * The template is compiled, argument binding happens and then the query is executed.
     */
    template <class... Types, class... ArgTypes>
    std::unique_ptr<CppSQLite3IteratorTable<Types...>>
    queryTable( const SQLQueryTemplateString<Types...>& rxQueryTemplate, const ArgTypes&... args )
    {
        return getTable<Types...>( rxQueryTemplate( *this ), args... );
    } // public method
}; // class

/* Creates the text for a SQL INSERT statement.
 * In the current version we automatically insert a NULL for the first column.
 * (We assume that the first column contains the primary key.)
 */
std::string
sCreateSQLInsertStatementText( const char* pcTableName, const unsigned int uiNumberOfArguemnts, bool = true );

/* Be warned: After the database object has been gone don't use the insert anymore.
 * Design flaw: ALthough the database has been gone the insert object can still survive. :-((
 */
class CppSQLiteExtStatementParent : public CppSQLite3Statement
{
  protected:
    /* Disallow copying an CppSQLiteExtStatementParent.
     * Copying of statements might have strange effects.
     */
    CppSQLiteExtStatementParent( const CppSQLiteExtStatementParent& ) = delete;

    /* Base-case of the binding of arguments. */
    template <int N> void bindArguments( )
    {} // base case

    /* Recursive case of the binding process */
    template <int N, typename First, typename... Rest> void bindArguments( const First& arg, const Rest&... args )
    {
        /* bind an argument using the corresponding method in CppSQLite3Statement
         */
        bind( N, arg );

        bindArguments<N + 1, Rest...>( args... );
    } // variadic method

    template <int N> void bindArgumentsForwarding( )
    {} // base case

    template <int N, typename First, typename... Rest> void bindArgumentsForwarding( First&& arg, Rest&&... args )
    {
        /* bind an argument using the corresponding method in CppSQLite3Statement
         */
        bind( N, std::forward<First>( arg ) );

        bindArgumentsForwarding<N + 1, Rest...>( std::forward<Rest>( args )... );
    } // variadic method

  public:
#if DEBUG_LEVEL > 0
    std::string sStatementText;
#endif

    /* ArgTypes are the types of the arguments of the query.
     * Return value is the number of rows changed by the statement.
     */
    template <typename... ArgTypes> int bindAndExecute( const ArgTypes&... args )
    {
        for( size_t uiTryCounter = 0; uiTryCounter < 2; uiTryCounter++ )
        {
            try
            { /* TO DO: Check whether this reset is really necessary. */
                reset( );

                /* This statement must be place here, check why. */
                bindArguments<1, ArgTypes...>( args... );

                /* Execute the single statement; this results in a fresh row in the corresponding table.*/
                return CppSQLite3Statement::execDML( );
            } // try
            catch( CppSQLite3Exception& xException )
            {
                /* Something went wrong, in most cases a timeout due to some lock.
                 * We try once again...
                 */
                std::cerr << "Database operation failed: [" << xException.errorMessage( ) << "] We try once again ..."
                          << std::endl;
                continue;
            } // catch
        } // for

        /* If you come to this point the operation failed 10 times.
         * This is critical and we give up.
         */
        throw CppSQLite3Exception( 1, "Database operation failed after 2 tries.", false );
    } // operator

    CppSQLiteExtStatementParent( CppSQLiteDBExtended& rxDatabase, const char* pcStatementText )
#if DEBUG_LEVEL > 0
        : sStatementText( pcStatementText )
#endif
    { /* We compile the statement by using the database object.
       * tricky use of the assignment operator defined in the base class. (hmm... C++)
       */
        CppSQLite3Statement::operator=( rxDatabase.compileStatement( pcStatementText ) );
    } // constructor

    /* Virtual destructor for inheritance purposes.
     */
    virtual ~CppSQLiteExtStatementParent( )
    {} // destructor
}; // class


/* SQL Query statement as specialization of a general SQL statement.
 * ...Types are the types of the columns of the returned table.
 */
template <class... Types> class CppSQLiteExtQueryStatementParent : public CppSQLiteExtStatementParent
{
  private:
    /* FIXME: It would be more secure to use a shared pointer over here. */
    CppSQLiteDBExtended& rxDatabase; // SHARED POINTER REQUIRED !!!!

  public:
    /* Constructor, where we get the SQL-statement from some query template. */
    // CppSQLiteExtQueryStatementParent<Types...>( CppSQLiteDBExtended& rxDatabase,
    //                                       const SQLQueryTemplateString<Types...>& rSQLQueryTemplate )
    //     : rxDatabase( rxDatabase ), CppSQLiteExtStatementParent( rxDatabase, rSQLQueryTemplate.pcText )
    // {} // constructor

    /* Constructor, where we get the SQL-statement directly as const string. */
    CppSQLiteExtQueryStatementParent<Types...>( CppSQLiteDBExtended& rxDatabase, const char* pcStatementText )
        : CppSQLiteExtStatementParent( rxDatabase, pcStatementText ), rxDatabase( rxDatabase )
    {} // constructor

    virtual ~CppSQLiteExtQueryStatementParent( )
    {} // destructor

    class Iterator
    {
      private:
        CppSQLite3Query xQuery;
        /* Repeatedly used functor object for column value extraction. */
        GetTupleElement xGetElementFunctor;

      public:
        template <class... ArgTypes>
        Iterator( CppSQLiteExtQueryStatementParent& rQuery, ArgTypes&&... args )
            : xQuery( rQuery.bindAndExecQuery( std::forward<ArgTypes>( args )... ) ), xGetElementFunctor( this->xQuery )
        {
            /* We check whether the number of columns is ok.
             * (Types are the types of the columns of the result table)
             */
            if( xQuery.numFields( ) != ( sizeof...( Types ) ) )
            {
                throw CppSQLite3Exception( CPPSQLITE_ERROR,
                                           "actual number of columns does not match number of types of template",
                                           false // => DONT_DELETE_MSG
                );
            } // if
        } // constructor

        std::tuple<Types...> get( )
        {
            std::tuple<Types...> SingleResultRowAsTuple;
            assert( !xQuery.eof( ) );
            iterateOverTupleCurry2( xGetElementFunctor, SingleResultRowAsTuple );
            return SingleResultRowAsTuple;
        } // method

        void next( )
        {
            xQuery.nextRow( );
        } // method

        bool eof( )
        {
            return xQuery.eof( );
        } // method
    }; // class

    /* ArgTypes are the types of the arguments of the query.
     * The class-types Types are the types of the columns of the returned table.
     * FIX ME: Use forwarding over here.
     */
    template <typename... ArgTypes> CppSQLite3Query bindAndExecQuery( const ArgTypes&... args )
    {
        for( size_t uiTryCounter = 0; uiTryCounter < 2; uiTryCounter++ )
        {
            try
            { /* TO DO: Check whether this reset is really necessary.
               * Reset the current statement.
               */
                reset( );

                /* This statement must be place here, check why. */
                bindArguments<1, ArgTypes...>( args... );

                /* Execute the single statement; this results in a fresh row in the corresponding table. */
                return CppSQLite3Statement::execQuery( );
            } // try
            catch( CppSQLite3Exception& xException )
            {
                /* Something went wrong, in most cases a timeout due to some lock.
                 * We try once again...
                 */
                std::cerr << "Database operation failed: [" << xException.errorMessage( ) << "] We try once again ..."
                          << std::endl;
                continue;
            } // catch
        } // for

        /* If you come to this point the operation failed 10 times.
         * This is critical and we give up.
         */
        throw CppSQLite3Exception( 1, "Database operation failed after 2 tries.", false );
    } // operator

    /* Creates a SQL-table (vector of tuples) as output.
     * Types shouldn't comprise some const char *, because of memory allocation problems
     */
    template <class... ArgTypes> std::unique_ptr<CppSQLite3IteratorTable<Types...>> exec( ArgTypes&&... args )
    {
        return rxDatabase.getTable( *this, std::forward<ArgTypes>( args )... );
    } // public method

    template <class... ArgTypes> Iterator vExecuteAndReturnIterator( ArgTypes&&... args )
    {
        return Iterator( *this, args... );
    } // method

    /* Applies the given function to all rows of the table.
     * May throw an exception if something goes wrong.
     */
    template <class TP_FUNC_APPLY, class... ArgTypes>
    void vExecuteAndForAllRowsDo(
        TP_FUNC_APPLY&& functor, // function called for each row retrived in the context of the query
        ArgTypes&&... args // query arguments
    )
    { /* We bind all arguments and execute the statement.
       * This process can fail and throw an exception.
       */
        CppSQLite3Query xQuery = this->bindAndExecQuery( std::forward<ArgTypes>( args )... );

        /* We check whether the number of columns is ok.
         * (Types are the types of the columns of the result table)
         */
        if( xQuery.numFields( ) != ( sizeof...( Types ) ) )
        {
            throw CppSQLite3Exception( CPPSQLITE_ERROR,
                                       "actual number of columns does not match number of types of template",
                                       false // => DONT_DELETE_MSG
            );
        } // if

        /* Repeatedly used functor object for column value extraction. */
        auto xGetElementFunctor = GetTupleElement( xQuery );

        /* Repeatedly used tuple temporary that keeps the content of a single row of the result table.
         */
        std::tuple<Types...> SingleResultRowAsTuple;

        /* We iterate over all rows of our table.
         */
        while( !xQuery.eof( ) )
        { /* Fill the tuple with the values of the columns in the current row. */
            iterateOverTupleCurry2( xGetElementFunctor, SingleResultRowAsTuple );

            /* Call function and deliver a parameter pack in form of a tuple. */
            functor( SingleResultRowAsTuple );

            /* Move to the next row of the table.
             */
            xQuery.nextRow( );
        } // while
    } // method

    /* Like vExecuteAndForAllRowsDo, but with the extension that the
     * tuple is unapcked to parameter by reference
     */
    template <typename TP_FUNC_APPLY, class... ArgTypes>
    void vExecuteAndForAllRowsUnpackedDo( TP_FUNC_APPLY&& functor, //
                                          ArgTypes&&... args )
    { /* We use vExecuteAndForAllRowsDo and unpack the tuples in the context of the functor calls.
       */
        vExecuteAndForAllRowsDo(
            [&functor]( const std::tuple<Types...>&
                            rxRowAsTuple ) { /* Unpack the tuple that contains the data of the current table-row and
                                              * call the functor given the unpacked tuple data as arguments.
                                              * Seems that perfect forwarding within a lambda works well, although there
                                              * is some capturing by reference involved.
                                              */
                                             ( TupleUnpackToParameterByReference<Types...>(
                                                 std::forward<TP_FUNC_APPLY>( functor ) ) )( rxRowAsTuple );
            }, // lambda

            std::forward<ArgTypes>( args )... // arguments of the query
        ); // function call
    } // method

    /* Maps one column of query into vector. */
    template <size_t N, class... ArgTypes>
    std::vector<typename std::tuple_element<N, std::tuple<Types...>>::type>
    executeAndStoreInVector( ArgTypes&&... args )
    { /* Returned vector that receives column of table
       */
        //// using type = typename std::tuple_element<N, std::tuple<Args...>>::type;

        std::vector<typename std::tuple_element<N, std::tuple<Types...>>::type> xReturnedVector;

        /* Iterate over all table rows
         */
        vExecuteAndForAllRowsDo(
            [&xReturnedVector]( const std::tuple<Types...>& rxRowAsTuple ) {
                xReturnedVector.emplace_back( std::get<N>( rxRowAsTuple ) );
            }, // lambda
            std::forward<ArgTypes>( args )... // arguments of the query
        ); // method call

        return xReturnedVector;
    } // template method

    /* Maps one column of query into vector. */
    template <class... ArgTypes>
    std::vector<typename std::tuple<Types...>> executeAndStoreAllInVector( ArgTypes&&... args )
    { /* Returned vector that receives column of table
       */
        //// using type = typename std::tuple_element<N, std::tuple<Args...>>::type;

        std::vector<typename std::tuple<Types...>> xReturnedVector;

        /* Iterate over all table rows
         */
        vExecuteAndForAllRowsDo(
            [&xReturnedVector]( const std::tuple<Types...>& rxRowAsTuple ) {
                xReturnedVector.push_back( rxRowAsTuple );
            }, // lambda
            std::forward<ArgTypes>( args )... // arguments of the query
        ); // method call

        return xReturnedVector;
    } // template method

    /* Maps one column of query into set. */
    template <size_t N, class... ArgTypes>
    std::set<typename std::tuple_element<N, std::tuple<Types...>>::type> executeAndStoreInSet( ArgTypes&&... args )
    {
        std::set<typename std::tuple_element<N, std::tuple<Types...>>::type> xReturnedSet;

        /* Iterate over all table rows
         */
        vExecuteAndForAllRowsDo(
            [&xReturnedSet]( const std::tuple<Types...>& rxRowAsTuple ) {
                xReturnedSet.emplace( std::get<N>( rxRowAsTuple ) );
            }, // lambda
            std::forward<ArgTypes>( args )... // arguments of the query
        ); // method call

        return xReturnedSet;
    } // template method

    /* Executes the query and applies fFunction to the first row, if this row does exist. */
    template <class FunctorType, class... ArgTypes>
    void vExecAndApplyToFirstRow( FunctorType&& fFunction, const ArgTypes&... args )
    {
        auto pTableRef = exec( args... );

        if( pTableRef->size( ) > 0 )
        {
            fFunction( *( pTableRef->begin( ) ) );
        } // if
    } // public method

    /* In the case of some scalar, we take the type of the first column as the type of the scalar.
     */
    typedef typename std::tuple_element<0, std::tuple<Types...>>::type ScalarType;

    /** Executes the current query, where the result is assumed to be a scalar.
     * Gets an function as argument, which is executed in the case that the query delivers a non-empty result.
     */
    template <class FunctorType, class... ArgTypes> void vScalar( FunctorType&& fFunction, ArgTypes&&... args )
    {
        auto pTableRef = exec( std::forward<ArgTypes>( args )... );

        if( pTableRef->size( ) > 0 )
        {
            fFunction( std::get<0>( *( pTableRef->begin( ) ) ) );
        } // if
    } // public method

    /** Executes the current query, where the result is assumed to be a scalar. Returns the scalar.
     * Check, whether there is some efficient way for implementing this method
     */
    template <class... ArgTypes> ScalarType scalar( ArgTypes&&... args )
    {
        auto pTableRef = exec( std::forward<ArgTypes>( args )... );

        if( pTableRef->size() == 0 )
            throw std::runtime_error( "EoF; no scalar value in query" );

        return std::get<0>( *( pTableRef->begin( ) ) );
    } // public method
}; // class


/**
 * @brief This class is no more than a wrapper around CppSQLiteExtStatementParent
 * @details
 * It is necessary, so that we can run EXPLAIN QUERY PLAN in debug mode.
 * Look at CppSQLiteExtStatementParent for all relevant method definitions & implementations
 */
template <class STATEMENT> class CppSQLiteExtDebugWrapper : public STATEMENT
{
  public:
#if DEBUG_LEVEL > 0
    /* In debug mode we create the ability to look at the output of the sqlite query planner:
     */
    CppSQLiteExtQueryStatementParent<int64_t, int64_t, int64_t, std::string> xExplainQueryPlan;
    bool bExplained = false;

    template <typename... ArgTypes> void bindAndExplain( const ArgTypes&... args )
    {
        if( bExplained )
            return;
        bExplained = true;

        std::cout << "Statement: " << STATEMENT::sStatementText << std::endl;
        std::cout << "id\tpId\t?\tdetail" << std::endl;
        std::cout << "--\t---\t-\t------" << std::endl;
        xExplainQueryPlan.vExecuteAndForAllRowsDo(
            []( auto xTuple ) {
                /* print result to std::cout */
                std::cout << std::get<0>( xTuple ) << "\t" << std::get<1>( xTuple ) << "\t" << std::get<2>( xTuple )
                          << "\t" << std::get<3>( xTuple ) << "\t" << std::endl;
            },
            args... );
        std::cout << std::endl;
    } // method
#endif

    CppSQLiteExtDebugWrapper( CppSQLiteDBExtended& rxDatabase, const char* pcStatementText )
        : STATEMENT( rxDatabase, pcStatementText )
#if DEBUG_LEVEL > 0
          ,
          xExplainQueryPlan( rxDatabase, std::string( "EXPLAIN QUERY PLAN " ).append( pcStatementText ).c_str( ) )
#endif
    {} // constructor

    /* Virtual destructor for inheritance purposes.
     */
    virtual ~CppSQLiteExtDebugWrapper( )
    {} // destructor
}; // class

class CppSQLiteExtStatement : public CppSQLiteExtDebugWrapper<CppSQLiteExtStatementParent>
{
  public:
    using CppSQLiteExtDebugWrapper<CppSQLiteExtStatementParent>::CppSQLiteExtDebugWrapper; // use constructor of
                                                                                           // superclass

#ifdef SQL_VERBOSE
#if DEBUG_LEVEL > 0
    /* @brief see CppSQLiteExtStatementParent::bindAndExecute
     */
    template <typename... ArgTypes> int bindAndExecute( const ArgTypes&... args )
    {
        bindAndExplain( args... );
        return CppSQLiteExtDebugWrapper<CppSQLiteExtStatementParent>::bindAndExecute( args... );
    } // method
#endif
#endif

    /* Virtual destructor for inheritance purposes.
     */
    virtual ~CppSQLiteExtStatement( )
    {} // destructor
}; // class

template <class... Types>
class CppSQLiteExtQueryStatement : public CppSQLiteExtDebugWrapper<CppSQLiteExtQueryStatementParent<Types...>>
{
  public:
    using CppSQLiteExtDebugWrapper<
        CppSQLiteExtQueryStatementParent<Types...>>::CppSQLiteExtDebugWrapper; // use constructor of superclass

#ifdef SQL_VERBOSE
#if DEBUG_LEVEL > 0
    /* @brief see CppSQLiteExtQueryStatementParent::bindAndExecute
     */
    template <typename... ArgTypes> int bindAndExecute( const ArgTypes&... args )
    {
        bindAndExplain( args... );
        return CppSQLiteExtDebugWrapper<CppSQLiteExtQueryStatementParent<Types...>>::bindAndExecute( args... );
    } // method
#endif
#endif

    /* Virtual destructor for inheritance purposes.
     */
    virtual ~CppSQLiteExtQueryStatement( )
    {} // destructor
}; // class

template <class... Types> struct SQLQueryTemplateString
{
    /* The request text itself.
     */
    const std::string sSQLQueryText;

    SQLQueryTemplateString( const char* const pcText ) : sSQLQueryText( pcText )
    {} // constructor

    /* Generator for a prepared statement stored in the database object.
     */
    CppSQLiteExtQueryStatement<Types...> operator( )( CppSQLiteDBExtended& rxDatabase ) const
    {
        return CppSQLiteExtQueryStatement<Types...>( rxDatabase, sSQLQueryText.c_str( ) );
    } // operator
}; // struct

/* Exception safe form of transactions for SQLite databases.
 */
class CppSQLiteExtImmediateTransactionContext
{
  private:
    /* The database for the transaction context
     */
    CppSQLiteDBExtended& rxDatabase;

  public:
    /* The constructor initializes the database and starts the transaction.
     */
    CppSQLiteExtImmediateTransactionContext( CppSQLiteDBExtended& rxDatabase ) : rxDatabase( rxDatabase )
    {
#if DEBUG_LEVEL > 0
        std::cout << "Begin transaction" << std::endl;
#endif
        rxDatabase.xStatementBeginTransaction->execDML( );
    } // constructor

    /* With the destruction of the object we end the transaction.
     */
    ~CppSQLiteExtImmediateTransactionContext( )
    {
#if DEBUG_LEVEL > 0
        std::cout << "End transaction" << std::endl;
#endif
        rxDatabase.xStatementEndTransaction->execDML( );
    } // destructor
}; // class CppSQLiteExtImmediateTransactionContext


template <typename TP_TYPE> std::string getSQLTypeName( )
{
    static_assert( std::is_base_of<SQL_BLOB, TP_TYPE>( ),
                   "Database column type must be either a primitive or inherit from SQL_BLOB!" );
    return "BLOB";
} // method

// list of valied types:
template <> std::string getSQLTypeName<std::string>( );
template <> std::string getSQLTypeName<bool>( );
// numeric:
// full number:
// sqLITE maps all numbers to INTEGER anyways?
// signed
template <> std::string getSQLTypeName<int8_t>( );
template <> std::string getSQLTypeName<int16_t>( );
template <> std::string getSQLTypeName<int32_t>( );
template <> std::string getSQLTypeName<int64_t>( );

// unsigned
template <> std::string getSQLTypeName<uint8_t>( );
template <> std::string getSQLTypeName<uint16_t>( );
template <> std::string getSQLTypeName<uint32_t>( );

// template <> std::string SQLTypeName<uint64_t>::name = "INTEGER"; <- this is too large for sqlite3

// floating point:
template <> std::string getSQLTypeName<double>( );
template <> std::string getSQLTypeName<float>( );


/* C++ wrapper/interface for a SQL-table.
 */
template <typename... Types> class CppSQLiteExtBasicTable
{
  protected:
    /* database that contains the table.
     * This construct can become quite dangerous. Better we take same shared pointer over here.
     */
    CppSQLiteDBExtended& rxDatabase;

    /* If this value is true we have an automatic primary key generation.
     */
    bool bAutomaticPrimaryKeyColumn;
    const bool bWithoutRowId;

    // rxDatabaseColumns = the names and type of the database columns.
    /* Create the table within the associated database.
     * May throw an exception if the process fails.
     */
    void vCreateTableInDB( const std::vector<std::string>& rxDatabaseColumns,
                           const std::vector<std::string>& vConstraints,
                           const bool bWithoutRowId )
    {
        /* Forward the creation of the table to the database object.
         * Improvement: Create a second statement, where we infer the column types in SQL from the C++ datatypes, given
         * in Types.
         */
        rxDatabase.vCreateTable(
            rxDatabaseColumns, sTableName.c_str( ), bAutomaticPrimaryKeyColumn, vConstraints, bWithoutRowId );
    } // method


    /**
     * @brief infer the type of each table column.
     */
    template <size_t N, typename TP_CURR, typename... TP_LIST> struct InferTypes
    {
        static void iterate( std::vector<std::string>& rxDatabaseColumns )
        {
            rxDatabaseColumns[ rxDatabaseColumns.size( ) - N ] += " " + getSQLTypeName<TP_CURR>( );
            InferTypes<N - 1, TP_LIST...>::iterate( rxDatabaseColumns );
        } // method
    }; // struct
    template <typename TP_CURR> struct InferTypes<1, TP_CURR>
    {
        static void iterate( std::vector<std::string>& rxDatabaseColumns )
        {
            rxDatabaseColumns.back( ) += " " + getSQLTypeName<TP_CURR>( );
        } // method
    }; // struct

  public:
    /* the name of the table
     */
    const std::string sTableName;

    /* Constructor that creates the table within the database
     */
    CppSQLiteExtBasicTable(
        CppSQLiteDBExtended& rxDatabase, // the database where the table resides
        const std::string& rsTableName, // name of the table in the database
        const std::vector<std::string>& rxDatabaseColumns, // column definitions of the table
        bool bAutomaticPrimaryKeyColumn, // true == we build automatically a column for the primary key
        const std::vector<std::string>& vConstraints = {},
        const bool bWithoutRowId = false )
        : rxDatabase( rxDatabase ),
          bAutomaticPrimaryKeyColumn( bAutomaticPrimaryKeyColumn ),
          bWithoutRowId( bWithoutRowId ),
          sTableName( rsTableName )
    {
        if( rxDatabase.eDatabaseOpeningMode == eCREATE_DB )
        {
            /* If the database is opened in a "creation mode", we create the table within the database.
             */
            std::vector<std::string> xTypedDatabaseColumns( rxDatabaseColumns );
            InferTypes<sizeof...( Types ), Types...>::iterate( xTypedDatabaseColumns );
            vCreateTableInDB( xTypedDatabaseColumns, vConstraints, bWithoutRowId );
        } // if
    } // constructor

    void clearTable( )
    {
        CppSQLiteExtQueryStatement<int64_t>( rxDatabase, std::string( "DELETE FROM " ).append( sTableName ).c_str( ) )
            .bindAndExecQuery<>( );
    } // method

    bool empty( )
    {
        int64_t iRes = CppSQLiteExtQueryStatement<int64_t>(
                           rxDatabase, std::string( "SELECT count(*) FROM " ).append( sTableName ).c_str( ) )
                           .scalar( );
        return iRes == 0;
    } // method
}; // class

/* SQL Insert statement as specialization of a general SQL statement.
 */
template <typename... Types> class CppSQLiteExtInsertStatement : public CppSQLiteExtStatement
{
    CppSQLiteDBExtended& rxDatabase;

  public:
    /* Constructor for insertion statement.
     */
    CppSQLiteExtInsertStatement( CppSQLiteDBExtended& rxDatabase, // the database where the insertion shall happen
                                 const char* pcTableName, // name of the table where we want to insert
                                 bool bFirstColumnAsNULL = true // shall the first column automatically inserted as NULL
                                                                // (for tables with automatic primary key)
                                 )
        : CppSQLiteExtStatement( rxDatabase,
                                 sCreateSQLInsertStatementText( pcTableName, sizeof...( Types ), bFirstColumnAsNULL )
                                     .c_str( ) ), // call superclass constructor
          rxDatabase( rxDatabase )
    {
        // std::cout << "Create Insertion statement for table \"" << pcTableName << "\"" << std::endl;
    } // constructor

    /* Destructor for insertion statement.
     */
    virtual ~CppSQLiteExtInsertStatement( )
    {} // destructor

    // void operator()( Types&& ... args )
    /**
     * returns the primary key of the inserted row... (if table has no primary key, the the rowid is returned)
     * FIX ME: Integrate some clean forwarding scheme.
     */
    int64_t operator( )( const Types&... args ) // deprecated
    {
        /* TO DO: Use meta-programming here
         */
        for( size_t uiTryCounter = 0; uiTryCounter < 10; uiTryCounter++ )
        {
            try
            {
                /* Bind all arguments for execution purposes
                 */
                // bindArgumentsForwarding<1, Types...>( std::forward<Types>( args )... );
                bindArguments<1, Types...>( args... ); // deprecated

                /* Execute the single statement; this results in a fresh row in the corresponding table.
                 */
                CppSQLite3Statement::execDML( );

                /* Reset the statement for the next execution.
                 * Reset was already done by execDML
                 */
                // reset();
            } // try
            catch( CppSQLite3Exception& xException )
            {
                if( xException.errorCode( ) == 19 ) // sqlite constraint failed
                    throw xException;
                /* Something went wrong, in most cases a timeout due to some lock.
                 * We try once again...
                 */
                std::cerr << "Database operation failed: [" << xException.errorMessage( ) << "] try once again ..."
                          << std::endl;
                continue;
            } // catch

            /* The insertion statement could be successfully executed ...
             * We indicate the end of the transaction.
             */
            return static_cast<int64_t>( rxDatabase.lastRowId( ) );
        } // for

        /* If you come to this point the operation failed 10 times.
         * This is critical and we give up.
         */
        throw CppSQLite3Exception( CPPSQLITE_ERROR,
                                   "Insertion failed finally after 10 tries",
                                   false // => DONT_DELETE_MSG
        );
    } // operator
}; // class

/* Generic function that simply serializes the arguments given and writes the serialized form to the output stream.
 */
template <typename Type> void genericPrintArgs( std::ostream& rxOutputStream, const Type& arg )
{
    /* bind an argument using the corresponding method in CppSQLite3Statement
     */
    rxOutputStream << arg << "\n";
} // variadic method (recursion end)

template <typename First, typename... Remaining>
void genericPrintArgs( std::ostream& rxOutputStream, const First& arg, const Remaining&... args )
{
    /* bind an argument using the corresponding method in CppSQLite3Statement
     */
    rxOutputStream << arg << " : ";

    /* Print the remaining arguments.
     */
    genericPrintArgs<Remaining...>( rxOutputStream, args... );
} // variadic method (recursive case)

/* WARNING: You have to declare a primary key, if you rely on this class for table representation.
 * Types are the types of table columns.
 * The types should match the types given to SQLite in the context of the table definition.
 */
template <typename... Types> class CppSQLiteExtTable : public CppSQLiteExtBasicTable<Types...>
{
  protected:
    /* In a parallel environment we have to serialize all database access.
     * This is done by the following mutex object, which is initialized as part of the constructor.
     */
    std::mutex xTableAccessMutex;

    /* The argument function is called for all rows of the table, where the column-values are delivered as actual
     * parameter values. If the argument is of type const char *&, the allocated memory of the string is available until
     * the end of the function.
     */
    template <typename... ArgTypes, typename F>
    void vForAllTableRowsUnpackedDoCore( F&& function
                                         //// ThreadPool *pxThreadPool = nullptr
    )
    { /* Create a temporary SELECT * FROM <table name> statement and execute it.
       * Then call function for each table row, with all columns as parameter of the function.
       */
        CppSQLiteExtQueryStatement<ArgTypes...>( this->rxDatabase,
                                                 std::string( "SELECT * FROM " ).append( this->sTableName ).c_str( ) )
            .vExecuteAndForAllRowsUnpackedDo( std::forward<F>( function ) /* , pxThreadPool */ );
    } // method

    /* Dumps the content of table alignments for debugging purposes to the given stream.
     * This core function can be used from subclasses using different type patterns.
     */
    template <typename... ArgTypes> void vDumpCore( std::ostream& rxOutputStream )
    {
        vForAllTableRowsUnpackedDoCore<ArgTypes...>(
            [&]( const ArgTypes&... args ) { genericPrintArgs<ArgTypes...>( rxOutputStream, args... ); } // lambda
        ); // function call
    } // method

    /* Buffer memory, that can be used for keeping table rows before insertion into the table.
     * Each table row is stored as single n-tuple in the vector.
     * Background: It boost performance tremendously to group several insertions into one single transaction.
     */
    std::vector<std::tuple<Types...>> xTableRowBuffer;

    /* Flushes the table buffer via a single transaction.
     * We expect here, that we manged already all task with respect to concurrent access.
     */
    void vFlushRowBufferCore( )
    { /* We start an immediate transaction. The transaction is automatically finished with the end of scope.
       * Doing so, we boost the performance, because all insertion statements can be executed in a single transaction.
       */
        //// std::cout << "BEGIN BUFFER FLUSH\n";
        CppSQLiteExtImmediateTransactionContext xBeginTransactionContext( this->rxDatabase );

        /* Do a buffer flush by inserting all buffered rows into the database.
         */
        for( std::tuple<Types...>& rxTableRow : xTableRowBuffer )
        {
            //// std::cout << rxTableRow << "\n";
            /* We apply the functor xInsertRow to the data of the unpacked tuple rxTableRow.
             * We have to avoid a copy of the functor xInsertRow or we are in big trouble.
             */
            TupleUnpackAndCallFunctor<Types...>( xInsertRow, rxTableRow );
        } // for
        //// std::cout << "END BUFFER FLUSH\n";
        /* end of method -> We automatically end the transaction context
         */
    } // method

  public:
    /* Statement object for inserting a fresh row into the table.
     * insertRow has to be called with argument types that match the definition of the table.
     * WARNING: Here is no mutex for synchronizing concurrent access.
     */
    CppSQLiteExtInsertStatement<Types...> xInsertRow;

    /* Stores a single table-row in the buffer for later flushing.
     * SAFE FOR CONCURRENT ACCESS TO TABLE!
     */
    void vInsertRowIntoBuffer( const Types&... args )
    { /* We have to use some mutex for serializing the access to the database if we are in concurrent environment.
       */
        std::lock_guard<std::mutex> lock( xTableAccessMutex );

        xTableRowBuffer.push_back( std::tuple<Types...>( args... ) );
    } // method

    /* Flushes the table-row buffer by inserting all rows into the table.
     * SAFE FOR CONCURRENT ACCESS TO TABLE!
     */
    void vFlushRowBuffer( )
    { /* We have to use some mutex for serializing the access to the database if we are in concurrent environment. */
        std::lock_guard<std::mutex> lock( xTableAccessMutex );

        vFlushRowBufferCore( );

        /* Erase all flushed rows in buffer.
         */
        xTableRowBuffer.clear( );
    } // method

    /* The argument function is called for all rows of the table, where the column-values are delivered as actual
     * parameter values. If the argument is of type const char *&, the allocated memory of the string is available until
     * the end of the function.
     */
    template <typename F> inline void vForAllTableRowsUnpackedDo( F&& functor )
    {
#if _MSC_VER
        /* Visual C++ 2013 crashes with the code of the else branch.
         */
        vForAllTableRowsUnpackedDoCore<Types...>( std::forward<F>( functor ) );
#else
        this->vForAllTableRowsUnpackedDoCore<Types...>( std::forward<F>( functor ) );
#endif
    } // method

    /* Dumps the content of table alignments for debugging purposes to the given stream. */
    inline void vDump( std::ostream& rxOutputStream )
    {
        this->vDumpCore<Types...>( rxOutputStream );
    } // method

    /* Table constructor. */
    CppSQLiteExtTable( CppSQLiteDBExtended& rxDatabase, // the database where the table resides
                       const std::string& rsTableName, // name of the table in the database
                       const std::vector<std::string>& rxDatabaseColumns, // column definitions of the table
                       bool bAutomaticPrimaryKeyColumn, // true == we build automatically a column for the primary key
                       const std::vector<std::string>& vConstraints = {},
                       const bool bWithoutRowId = false )
        : CppSQLiteExtBasicTable<Types...>( rxDatabase,
                                            rsTableName,
                                            rxDatabaseColumns,
                                            bAutomaticPrimaryKeyColumn,
                                            vConstraints,
                                            bWithoutRowId ), // call superclass constructor
          xTableAccessMutex( ), // initialize the mutex for synchronizing table access in a concurrent environment
          xTableRowBuffer( ), // initialize the table row buffer
          xInsertRow( rxDatabase,
                      CppSQLiteExtBasicTable<Types...>::sTableName.c_str( ),
                      CppSQLiteExtBasicTable<Types...>::bAutomaticPrimaryKeyColumn ) // precompile insert statement
    {} // constructor

    /* Table destructor.
     */
    ~CppSQLiteExtTable( )
    { /* Flush the row buffer, if not done yet. */
        try
        {
            vFlushRowBuffer( );
        } // try
        catch( std::exception& rxException )
        {
            std::cerr << "~CppSQLiteExtTable() failed due to : " << rxException.what( );
        } // catch
        catch( ... )
        {
            std::cerr << "~CppSQLiteExtTable() failed due to unknown reason";
        } // catch
    } // method
}; // class

/* WARNING: You shouldn't declare a primary key if you rely on this class for table representation.
 * Types are the types of table columns.
 * The types should match the types given to SQLite in the context of the table definition.
 */
template <typename... Types> class CppSQLiteExtTableWithAutomaticPrimaryKey : public CppSQLiteExtTable<Types...>
{
  public:
    /* The argument function is called for all rows of the table, where the column-values are delivered as actual
     * parameter values. If the argument is of type const char *&, the allocated memory of the string is available until
     * the end of the function.
     */
    void vForAllTableRowsUnpackedDo( const std::function<void( int, Types... )>& function )
    {
        /* See http://stackoverflow.com/questions/5802077/strange-gcc-error-expected-primary-expression-before-token
         */
        this->template vForAllTableRowsUnpackedDoCore<int, Types...>( function );
    } // method

    /* Dumps the content of table alignments for debugging purposes.
     * Dumps the content to a given stream.
     */
    void vDump( std::ostream& rxOutputStream )
    {
        /* See http://stackoverflow.com/questions/5802077/strange-gcc-error-expected-primary-expression-before-token
         */
#if _MSC_VER
        /* Visual C++ 2013 crashes with the below code
         */
        vDumpCore<int, Types...>( rxOutputStream );
#else
        this->template vDumpCore<int, Types...>( rxOutputStream );
#endif
    } // method

    /* Constructor for table with automatic primary key.
     */
    CppSQLiteExtTableWithAutomaticPrimaryKey(
        CppSQLiteDBExtended& rxDatabase, // the database where the table resides
        const std::string& rsTableName, // name of the table in the database
        const std::vector<std::string>& rxDatabaseColumns, // column definitions of the table
        const std::vector<std::string>& vConstraints = {} )
        :
#if _MSC_VER
          /* Visual C++ 2013 crashes with the below code
           */
          CppSQLiteExtTable( rxDatabase, rsTableName, rxDatabaseColumns, true, vConstraints, false )
#else
          CppSQLiteExtTable<Types...>( rxDatabase, rsTableName, rxDatabaseColumns, true, vConstraints, false )
#endif
    {} // constructor
}; // class

/* WARNING: WITHOUT ROWID optimization is only helpfull in very few scenarios
 * read https://www.sqlite.org/withoutrowid.html before using this tabletype
 */
template <typename... Types> class CppSQLiteExtTableWithoutRowId : public CppSQLiteExtTable<Types...>
{
  public:
    /* The argument function is called for all rows of the table, where the column-values are delivered as actual
     * parameter values. If the argument is of type const char *&, the allocated memory of the string is available until
     * the end of the function.
     */
    void vForAllTableRowsUnpackedDo( const std::function<void( Types... )>& function )
    {
        /* See http://stackoverflow.com/questions/5802077/strange-gcc-error-expected-primary-expression-before-token
         */
        this->template vForAllTableRowsUnpackedDoCore<Types...>( function );
    } // method

    /* Dumps the content of table alignments for debugging purposes.
     * Dumps the content to a given stream.
     */
    void vDump( std::ostream& rxOutputStream )
    {
        /* See http://stackoverflow.com/questions/5802077/strange-gcc-error-expected-primary-expression-before-token
         */
#if _MSC_VER
        /* Visual C++ 2013 crashes with the below code
         */
        vDumpCore<Types...>( rxOutputStream );
#else
        this->template vDumpCore<Types...>( rxOutputStream );
#endif
    } // method

    /* Constructor for table without rowid.
     */
    CppSQLiteExtTableWithoutRowId( CppSQLiteDBExtended& rxDatabase, // the database where the table resides
                                   const std::string& rsTableName, // name of the table in the database
                                   const std::vector<std::string>& rxDatabaseColumns, // column definitions of the table
                                   const std::vector<std::string>& vConstraints = {} )
        :
#if _MSC_VER
          /* Visual C++ 2013 crashes with the below code
           */
          CppSQLiteExtTable( rxDatabase, rsTableName, rxDatabaseColumns, false, vConstraints, true )
#else
          CppSQLiteExtTable<Types...>( rxDatabase, rsTableName, rxDatabaseColumns, false, vConstraints, true )
#endif
    {} // constructor
}; // class
#endif