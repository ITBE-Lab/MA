/* Authors: Arne Kutzner and Markus Schmidt
 * Created: July 2020
 * This file is part of the ITBE-Lab code collection.
 * MIT License
 * @postgre_sql_core.h
 * @brief Core part of the PostgreSQL engine support. This file is intended to be included via postgre_sql.h merely.
 */

// Required for MSVC for avoiding conflicts with std::min, std::max
#ifdef _MSC_VER
#define NOMINMAX
#endif

// Inform later imports about the SQL backend
#define POSTGRESQL

// Deliver comprehensive error messages
#define VERBOSE_ERROR_MESG

// STL imports
#include <atomic> // atomic fetch_add
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "db_common.h"

// Include PostgreSQL header
#include "libpq-fe.h"

// For getting the code working with gcc 6.x.x compiler
#if( __GNUC__ && ( __GNUC__ < 7 ) && !defined(__clang__) )
#include <experimental/tuple>
#define STD_APPLY std::experimental::apply
// Part of the C++17 standard now
// std::ostream& operator<<( std::ostream& s, std::nullptr_t )
// {
//     return s << "nullptr";
// } // operator overload
#else
#include <tuple>
#define STD_APPLY std::apply
#endif

#if (defined( __GNUC__ ) && ( __GNUC__ < 8 )) && !(defined(__clang__) && (__clang_major__ >= 10))
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

// Integrate json support
// MySQL Database configuration is described via JSON
#include <json.hpp>
using nlohmann::json;

// Constants for definitions via json
const std::string CONNECTION = "CONNECTION";
const std::string HOSTNAME = "HOSTNAME";
const std::string USER = "USER";
const std::string PASSWORD = "PASSWORD";
const std::string PORT = "PORT";
const std::string MYSQL_EXTRA = "MYSQL_EXTRA";
const std::string NO_RESULT_STORAGE_ON_CLIENT = "NO_RESULT_STORAGE_ON_CLIENT";
const std::string FLAGS = "FLAGS";

/* Example for definition of database connection in json:
 * auto jDBConfig = json{ { SCHEMA, { { NAME, "used_schema" }, { FLAGS, { DROP_ON_CLOSURE } } } },
 *                        { CONNECTION,
 *                          { { HOSTNAME, "localhost" },
 *                            { USER, "root" }, //
 *                            { PASSWORD, "admin" },
 *                            { PORT, 0 } } } };
 */
#if 0
// MySQL import
#include <mysql.h>
// The header <my_global.h> has been gone since MySQL 8.
// The below typedef keeps the code compiling for older releases of MySQL.
#if !defined( MARIADB_BASE_VERSION ) && !defined( MARIADB_VERSION_ID ) && MYSQL_VERSION_ID >= 80001 &&                 \
    MYSQL_VERSION_ID != 80002
typedef bool my_bool;
#endif
#endif

#define DEFAULT_HOSTNAME "localhost"
#define DEFAULT_USER "postgres"
#define DEFAULT_PASSWD "admin"
#define DEFAULT_PORT 0

#define LAMBDA_WITH_PARAMETER_PACK_SUPPORTED

/* If we're not using GNU C, elide __attribute__ */
#ifndef __GNUC__
#define __attribute__( x ) /* NOTHING */
#endif

/* Type correspondences:
 * C++              SQL            MySQL C-API
 * std::string    - LONGTEXT     - MYSQL_TYPE_LONG_BLOB
 * GenericBlob    - LONGBLOB     - MYSQL_TYPE_LONG_BLOB
 * bool		      - TINYINT      - MYSQL_TYPE_TINY
 * uint32_t		  - INT UNSIGNED -
 *                               - MYSQL_TYPE_GEOMETRY
 * std::nullptr_t - NULL         - MYSQL_TYPE_NULL
 */

// Required for workaround with respect to explicit specialization of a template function in a template class.
// See: https://stackoverflow.com/questions/3052579/explicit-specialization-in-non-namespace-scope
template <typename T> struct identity
{
    typedef T type;
}; // struct

// Corresponding C++ type for PG type BIGSERIAL
struct PGBigSerial
{
    uint64_t id;

    PGBigSerial( uint64_t id ) : id( id )
    {}
}; // struct

// typedef uint64_t PGBigSerial;


/** @brief Exception that informs about a problem with PostgreSQL
 */
class PostgreSQLError : public std::runtime_error
{
  private:
    /** @brief Generate a MySQL error message
     */
    std::string makeMessage( PGconn* pDBConn )
    {
        const char* pErrorText = pDBConn ? PQerrorMessage( pDBConn ) : NULL;

        std::stringstream xBuffer;
        xBuffer << "PostgreSQL DB error: " << ( pErrorText ? pErrorText : "UNKNOWN ERROR" ) << std::endl;
        return xBuffer.str( );
    } // method


  public:
    /** @brief MySQLConException exception
     */
    PostgreSQLError( const std::string& rsErrExplain ) : std::runtime_error( "PostgreSQL error: " + rsErrExplain )
    {}

    /** @brief  Ask MySQL for the error text
     */
    PostgreSQLError( PGconn* pDBConn ) : std::runtime_error( makeMessage( pDBConn ) )
    {}
}; // class

/** @brief Minimalistic generic Blob (pointer plus size)
 */
struct GenericBlob
{
    char* pBuf; // pointerToBuffer
    size_t uiSize; // size of buffer
}; // struct

struct PG_GLobalEnv
{
    std::unique_ptr<std::atomic<size_t>> uiStmtCounter = std::make_unique<std::atomic<size_t>>( 0 );

    size_t getStmtId( ) const
    {
        return uiStmtCounter->fetch_add( 1 );
    } // method
}; // struct

/** @brief Single master object for concurrency synchronization. */
extern PG_GLobalEnv xPG_GLobalEnv;

/** @brief Big endian to little endian translation for 8 bytes */
template <typename TYPE> TYPE byteswap8( char* pVal )
{
    union
    {
        TYPE tpVal;
        uint64_t uiVal;
    } buf;

#if( defined _MSC_VER )
    buf.uiVal = _byteswap_uint64( *( (uint64_t*)pVal ) ); // MSVC
#elif( defined __GNUC__ )
    buf.uiVal = __builtin_bswap64( *( (uint64_t*)pVal ) ); // GCC
#else
#error Please insert an efficent 8 byte swap for your compiler
#endif
    return buf.tpVal;
} // method

/** @brief Big endian to little endian translation for 4 bytes */
template <typename TYPE> TYPE byteswap4( char* pVal )
{
    union
    {
        TYPE tpVal;
        uint32_t uiVal;
    } buf;

#if( defined _MSC_VER )
    buf.uiVal = _byteswap_ulong( *( (uint32_t*)pVal ) ); // MSVC
#elif( defined __GNUC__ )
    buf.uiVal = __builtin_bswap32( *( (uint32_t*)pVal ) ); // GCC
#else
#error Please insert an efficent 8 byte swap for your compiler
#endif
    return buf.tpVal;
} // method

/** @brief Big endian to little endian translation for 2 bytes */
template <typename TYPE> TYPE byteswap2( char* pVal )
{
    union
    {
        TYPE tpVal;
        uint16_t uiVal;
    } buf;
    buf.uiVal = *( (uint16_t*)pVal );
    buf.uiVal = ( buf.uiVal >> 8 ) | ( buf.uiVal << 8 );
    return buf.tpVal;
} // method

/* Basic class for a single cell in row for a query outcome.
 * Full support of C++17 would allow moving all these specializations inside MySQLConDB.
 */
template <typename CellType> class PGRowCellBase
{
  public:
    CellType* pCellValue; // pointer to the actual cell value (this pointer refers into tCellValues)
    size_t uiColNum; // column number of cell in query
    bool isNull = false; // if true, cell keeps a null value

    /** @brief Initialization of a cell must be done via this init method. */
    inline void init( CellType* pCellValue, // pointer to actual cell value
                      size_t uiColNum ) // column number of cell in query outcome
    {
        this->pCellValue = pCellValue; // backup pointer to cell value
        this->uiColNum = uiColNum; // backup column number
    } // method

    /** @brief Gets the value pointer for the current cell */
    inline char* getValPtr( const PGresult* pPGRes )
    {
        return PQgetvalue( pPGRes, 0, (int)this->uiColNum );
    } // method

    /** @brief Gets the length of the current cell in the case of binary data */
    inline int getValLength( const PGresult* pPGRes )
    {
        return PQgetlength( pPGRes, 0, (int)this->uiColNum );
    } // method

    // Dummy method that should be redefined by derived classes.
    inline void store( const PGresult* pPGRes );
    // {
    //     throw PostgreSQLError( std::string( "store in PGRowCellBase called for type: " ) + typeid( CellType ).name( )
    //     );
    // } // method
}; // class (PGRowCellBase)

/* This class is intended to be specialized for specific types. */
template <typename CellType> struct PGRowCell : public PGRowCellBase<CellType>
{}; // class

/* The trick here is the overriding of functions in the specialized subclasses.
 * Important: This explicit specializations must occur before any possible implicit specializations or
 *            the compiler complains about "re-specialization".
 */

// boolean value (In PG stored via 16 bit integers)
template <> struct PGRowCell<bool> : public PGRowCellBase<bool>
{
    inline void init( bool* pCellValue, size_t uiColNum )
    {
        PGRowCellBase::init( pCellValue, uiColNum );
    } // method

    inline void store( const PGresult* pPGRes )
    {
        *pCellValue = byteswap2<bool>( this->getValPtr( pPGRes ) );
        this->isNull = false;
    } // method
}; // specialized class

// 32 bit integer
template <> struct PGRowCell<int32_t> : public PGRowCellBase<int32_t>
{
    inline void init( int32_t* pCellValue, size_t uiColNum )
    {
        PGRowCellBase::init( pCellValue, uiColNum );
    } // method

    inline void store( const PGresult* pPGRes )
    {
        *pCellValue = byteswap4<int32_t>( this->getValPtr( pPGRes ) );
        this->isNull = false;
    } // method
}; // specialized class

// unsigned 32 bit integer are represented via 64 bit integers, because PG does not know unsigned types.
template <> struct PGRowCell<uint32_t> : public PGRowCellBase<uint32_t>
{
    inline void init( uint32_t* pCellValue, size_t uiColNum )
    {
        PGRowCellBase::init( pCellValue, uiColNum );
    } // method

    inline void store( const PGresult* pPGRes )
    {
        *pCellValue = (uint32_t)byteswap8<int64_t>( this->getValPtr( pPGRes ) );
        this->isNull = false;
    } // method
}; // specialized class

// 64 bit integer
template <> struct PGRowCell<int64_t> : public PGRowCellBase<int64_t>
{
    inline void init( int64_t* pCellValue, size_t uiColNum )
    {
        PGRowCellBase::init( pCellValue, uiColNum );
    } // method

    inline void store( const PGresult* pPGRes )
    {
        *pCellValue = byteswap8<int64_t>( this->getValPtr( pPGRes ) );
        this->isNull = false;
    } // method
}; // specialized class

//  unsigned 64 bit integer (unsigned info is extracted inside PGRowCellBase::init )
template <> struct PGRowCell<uint64_t> : public PGRowCellBase<uint64_t>
{
    inline void init( uint64_t* pCellValue, size_t uiColNum )
    {
        PGRowCellBase::init( pCellValue, uiColNum );
    } // method

    inline void store( const PGresult* pPGRes )
    {
        *pCellValue = (uint64_t)byteswap8<int64_t>( this->getValPtr( pPGRes ) );
        this->isNull = false;
    } // method
}; // specialized class

// 64 bit double
template <> struct PGRowCell<double> : public PGRowCellBase<double>
{
    inline void init( double* pCellValue, size_t uiColNum )
    {
        PGRowCellBase::init( pCellValue, uiColNum );
    } // method

    inline void store( const PGresult* pPGRes )
    {
        *pCellValue = byteswap8<double>( this->getValPtr( pPGRes ) );
        this->isNull = false;
    } // method
}; // specialized class

// std::string
template <> struct PGRowCell<std::string> : public PGRowCellBase<std::string>
{
    inline void init( std::string* pCellValue, size_t uiColNum )
    {
        PGRowCellBase::init( pCellValue, uiColNum );
    } // method

    inline void store( const PGresult* pPGRes )
    {
        // The binary representation of TEXT is, well, text, and since libpq was nice enough to append a zero byte to
        // it, it'll work just fine as a C string.
        pCellValue->assign( this->getValPtr( pPGRes ) );
    } // method
}; // specialized class

#if 1
/** @brief Models an "immediate" (not prepared) query. This kind of queries is fully parsed with each execution.
 */
template <typename DBConPtrType> class ImmediateQueryTmpl
{
  private:
    DBConPtrType pDBConn; // keeps MySQL connection
    std::string sQueryTxt; // statement text of query
    PGresult* pRes; // Query result

    /** @brief Make PostgreSQL error as text. */
    std::string errMsg( )
    {
        std::string sErrTxt = "Immediate query execution failed:\n";
#if 1
        sErrTxt.append( "Affected query: " ).append( sQueryTxt ).append( "\n" );
        sErrTxt.append( "Detailed PostgreSQL error msg:\n" );
#endif
        sErrTxt.append( PQerrorMessage( pDBConn->pPGConn ) ).append( "\n" );
        return sErrTxt;
    } // method

    /** @brief Must be called for avoiding memory leaks. */
    void clearResult( )
    {
        if( this->pRes ) // avoid memory leak
        {
            PQclear( this->pRes );
            this->pRes = NULL;
        } // if
    } // method

  public:
    /** @brief Constructs an immediate query. */
    ImmediateQueryTmpl( DBConPtrType pDBConn, const std::string& rsQueryTxt )
        : pDBConn( pDBConn ), sQueryTxt( rsQueryTxt ), pRes( NULL )
    {}

    /** @brief Executes the query and stores the result on the client. */
    SimpleTextTable exec( )
    {
        clearResult( );
        // Execute query
        this->pRes = PQexec( pDBConn->pPGConn, this->sQueryTxt.c_str( ) );
        if( PQresultStatus( this->pRes ) != PGRES_TUPLES_OK )
        {
            clearResult( );
            throw PostgreSQLError( this->errMsg( ) );
        } // if

        // Get all column names
        std::vector<std::string> vFieldNames;
        int iNumFields = PQnfields( this->pRes );
        for( int i = 0; i < iNumFields; i++ )
        {
            vFieldNames.emplace_back( PQfname( this->pRes, i ) );
        }

        // Get the rows of the returned table
        SimpleTextTable xretTbl( vFieldNames );
        for( int i = 0; i < PQntuples( this->pRes ); i++ )
        {
            std::vector<std::string> vRow;
            for( int j = 0; j < iNumFields; j++ )
                vRow.emplace_back( PQgetvalue( this->pRes, i, j ) );
            xretTbl.addRow( vRow );
        } // for
        clearResult( );

        return xretTbl;
    } // method

    /** @brief Frees all allocated resources. */
    ~ImmediateQueryTmpl( )
    {
        if( this->pRes )
            PQclear( this->pRes );
    } // destructor
}; // class ( ImmediateQueryTmpl )
#endif

/** @brief Explains MySQL queries and statements. Required for query optimization purposes. */
template <typename DBConPtrType> class QueryExplainer
{
    DBConPtrType pDBConn; // keeps MySQL connection

  public:
    /** @brief Explains the execution strategy for a MySQL statement.
     *  The statement delivered for inspection must be free of placeholders.
     */
    void explain( std::string rsStmtText )
    {
        // Create and execute EXPLAIN-statement.
        std::string sExplStmtText = "EXPLAIN " + rsStmtText;
        ImmediateQueryTmpl<DBConPtrType> xExplainQuery( pDBConn, sExplStmtText );
        auto xTbl = xExplainQuery.exec( );
        xTbl.print( );
    } // method

    /** @brief Constructs query explainer using the MySQL connection given as argument */
    QueryExplainer( DBConPtrType pDBConn ) : pDBConn( pDBConn )
    {}
}; // class ( QueryExplainer )


/** @brief MqSQL connection to database
 *  Warning: Class is not thread-safe.
 */
class PostgreSQLDBCon
{
  public:
    /// @brief Translates C++ to PostgresSQ Column Types
    class TypeTranslator
    {
      public:
        /* Translation of C++ type to SQLite column types.
         * Uses the specialization to overloading trick here.
         */
        template <typename Type> static inline std::string getSQLTypeName( )
        {
            return getSQLTypeName( identity<Type>( ) );
        } // method

        /** @brief Delivers the appropriate string for placeholder in prepared statement (as e.g. INSERT statement).
         * In most cases this is simply a '?'. But for some types, as e.g. spatial types, an additional wrapping is
         * required for indicating the precise kind of data-representation passed to the MySQL client. In such cases an
         * additional specialization is used for performing the wrapping.
         */
        template <typename Type> static inline std::string getPlaceholderForType( const std::string& rsInsertedText )
        {
            return rsInsertedText;
        } // method

      private:
        // list of valid types:
        // https://www.postgresql.org/docs/current/datatype.html
        template <typename Type> static inline std::string getSQLTypeName( identity<Type> );

        static inline std::string getSQLTypeName( identity<PGBigSerial> )
        {
            return "bigserial"; // used for primary key definitions merely
        } // private method

        // https://www.postgresql.org/docs/current/datatype-character.html
        static inline std::string getSQLTypeName( identity<std::string> )
        {
            return "text";
        } // private method

        // https://www.postgresql.org/docs/current/datatype-binary.html#DATATYPE-BINARY-TABLE
        static inline std::string getSQLTypeName( identity<GenericBlob> )
        {
            return "bytea";
        } // private method

        static inline std::string getSQLTypeName( identity<bool> )
        {
            return "int2";
        } // private method

        static inline std::string getSQLTypeName( identity<int32_t> )
        {
            return "int4";
        } // private method

        // WARNING: No support of unsigned types in PostgreSQL
        static inline std::string getSQLTypeName( identity<uint32_t> )
        {
            return "int8";
        } // private method

        static inline std::string getSQLTypeName( identity<int64_t> )
        {
            return "int8";
        } // private method

        // WARNING: No support of unsigned types in PostgreSQL
        static inline std::string getSQLTypeName( identity<uint64_t> )
        {
            return "int8";
        } // private method

        static inline std::string getSQLTypeName( identity<double> )
        {
            return "float8";
        } // private method

        static inline std::string getSQLTypeName( identity<nullptr_t> )
        {
            throw PostgreSQLError( "PG: Invalid request in getSQLTypeName for the type nullptr_t" );
            return "null";
        } // private method
    }; // class


#define PG_TEXT_ARG 0
#define PG_BINARY_ARG 1

    /** @brief Input argument of a prepared MySQL query. */
    class StmtArg
    {
      private:
        // Specifies the actual value of the parameters. A null pointer means the corresponding parameter
        // is null; otherwise the pointer points to a zero-terminated text string (for text format) or binary data in
        // the format expected by the server (for binary format).
        char*& rpParamValue;

        // Specifies the actual data length of a binary-format parameter. It is ignored for null parameters and
        // text-format parameters.
        int& riParamLength;

        // Specifies whether parameters are text (put a zero in the array entry for the corresponding
        // parameter) or binary (put a one in the array entry for the corresponding parameter).
        int& riParamFormat;

        std::string sTxtBuf; // used for textual representation of numbers

        // translation of number to string without using scientific notation
        template <typename T> inline std::string toString( T val )
        {
            std::stringstream xOutStream;
            xOutStream << std::fixed;
            xOutStream << val;
            return xOutStream.str( );
        } // method

        // const std::map<std::string, Oid>& rmOidMap;

        // We assume that all calls use the same type and cache this type for avoiding repeated map lookups
        // int getOid( const std::string& rTypeName )
        // {
        //     if( riParamFormat == std::numeric_limits<int>::min( ) )
        //         return riParamFormat = rmOidMap.at(rTypeName);
        //     return riParamFormat;
        // } // method

      public:
        // buffer for length in MYSQL_BIND (32 Bit on MSVC, 64 Bit on GCC)
        // Only required for BLOB and LONGTEXT
        unsigned long uiLength;

        StmtArg( char*& rpParamValue, int& riParamLength, int& riParamFormat )
            : rpParamValue( rpParamValue ), riParamLength( riParamLength ), riParamFormat( riParamFormat )
        {
            // By default we assume a text parameter
            riParamLength = 0;
            riParamFormat = PG_TEXT_ARG;
        } // constructor

        /* templated dispatcher.
         * The dispatcher allows using template specializations as replacement for overloading.
         * TIP: This dispatcher can be specialized as well. (See e.g. for NucSeq)
         */
        template <typename Type> inline void set( const Type& dArg )
        {
            // DEBUG: std::cout << "set called with dArg type: " << typeid( dArg ).name( ) << std::endl;
            setOverloaded( dArg );
        } // method

      private:
        /** @brief Binder for arguments in the context of statement execution:
         *  Remainder: You have to guarantee the lifetime of the data bound!
         */
        void inline setOverloaded( void* const& p ) = delete; // pointers are not allowed yet

        void inline setOverloaded( const bool& bVal )
        {
            sTxtBuf = toString( int( bVal ) );
            rpParamValue = (char*)sTxtBuf.c_str( );
        } // method

        void inline setOverloaded( const int32_t& iVal )
        {
            sTxtBuf = toString( iVal );
            rpParamValue = (char*)sTxtBuf.c_str( );
        } // method

        void inline setOverloaded( const uint32_t& iVal )
        {
            sTxtBuf = toString( iVal );
            rpParamValue = (char*)sTxtBuf.c_str( );
        } // method

        void inline setOverloaded( const int64_t& iVal )
        {
            sTxtBuf = toString( iVal );
            rpParamValue = (char*)sTxtBuf.c_str( );
        } // method

        void inline setOverloaded( const uint64_t& iVal )
        {
            sTxtBuf = toString( iVal );
            rpParamValue = (char*)sTxtBuf.c_str( );
        } // method

        void inline setOverloaded( const double& dVal )
        {
            sTxtBuf = toString( dVal );
            rpParamValue = (char*)sTxtBuf.c_str( );
        } // method

        // Corresponding head for C-NULL pointer:
        //
        void inline setOverloaded( const std::nullptr_t& p ) // NULL
        {
            rpParamValue = (char*)NULL;
            riParamLength = 0;
            riParamFormat = PG_BINARY_ARG;
        } // method

        // This should work on most architectures.
        // Due to the C++ standard each character of a string is a single byte.
        void inline setOverloaded( const std::string& rsText ) // Text
        {
            rpParamValue = (char*)rsText.c_str( );
            riParamLength = 0;
            riParamFormat = PG_TEXT_ARG;
        } // method

        void inline setOverloaded( const GenericBlob& rxBlob ) // Binary Blob
        {
            rpParamValue = rxBlob.pBuf;
            riParamLength = (int)rxBlob.uiSize;
            riParamFormat = PG_BINARY_ARG;
        } // method
    }; // inner class StmtArg

    /** @brief Prepared MySQL statement
     */
    template <typename DBPtrType> class PreparedStmtTmpl
    {
      protected:
        std::string sStmtName; // identifies the statement and has to be set by the constructor
        std::vector<char*> vParamValues; // vector of pointers keeping parameter value
        std::vector<int> vparamLengths; // vector of ints keeping parameter lengths
        std::vector<int> vparamFormats; // vector of ints keeping parameter formats
        int iStmtNumParams; // number of parameters of statement
        std::vector<StmtArg> vInputArgs; // the array of wrapped binders with references to the MySQL Bind
        DBPtrType const pDBConn; // DB Connection that the statement is embedded in
        PGresult* pPGRes = NULL; // Result of PG stmt execution
        bool bUseAsyncQuery; // If true, use PQsendQueryPrepared instead of PQexecPrepared
        const std::string sStmtText; // keep the statement text for more comprehensive error messages.

        /** @brief Translates a C++ type to its Oid via the corresponding PG SQL type
         */
        template <typename Type> Oid inline getOidForType( const Type& bVal )
        {
            return this->pDBConn->mOidMap.at( TypeTranslator::template getSQLTypeName<Type>( ) );
        } // method

        // Base case for binding with forwarding
        template <int N> void bindArgumentsForwarding( )
        {} // base case
        // Recursive case for binding with forwarding
        template <int N, typename Type, typename... Rest> void bindArgumentsForwarding( Type&& arg, Rest&&... args )
        {
            // bind argument N
            this->vInputArgs[ N ].set( std::forward<Type>( arg ) );
            // recursively bind remaining arguments
            bindArgumentsForwarding<N + 1, Rest...>( std::forward<Rest>( args )... );
        } // variadic method

        // Base case for binding with forwarding
        template <typename... Rest> void bindArgumentsForwardingRuntime( size_t N )
        {} // base case

        /** @brief as bindArgumentsForwarding, but using a runtime offset */
        template <typename Type, typename... Rest>
        void bindArgumentsForwardingRuntime( size_t N, Type&& arg, Rest&&... args )
        {
            // bind argument N
            this->vInputArgs[ N ].set( std::forward<Type>( arg ) );
            // recursively bind remaining arguments
            bindArgumentsForwardingRuntime<Rest...>( N + 1, std::forward<Rest>( args )... );
        } // variadic method

        /** @brief Returns a string that contains info regarding types and values of an argument pack.
         *  Primarily for debugging purposes...
         */
        template <typename... ArgTypes> void argsOIDs( ArgTypes&&... args )
        {
            std::vector<Oid> vOids;
            ( vOids.push_back( getOidForType( args ) ), ... ); // C++17

            // DEBUG: for( auto i : vOids )
            // DEBUG:     std::cout << "arg:" << i << std::endl;
        } // variadic method

        /** @brief Checks status of result stmt and throw exception in case of failure */
        void checkResultStatus( )
        {
            // FIXME: Isn't here a memory leak,because we do not clear pPGRes
            if( ( PQresultStatus( this->pPGRes ) == PGRES_COMMAND_OK ||
                  PQresultStatus( this->pPGRes ) == PGRES_TUPLES_OK ) )
                return;
            throwPostgreSQLError( );
        } // method

        /** @brief Throws an PostgreSQLError */
        void throwPostgreSQLError( )
        {
            throw PostgreSQLError( std::string( "Status Code:" ) + std::to_string( PQresultStatus( this->pPGRes ) ) +
                                   "\n" + std::string( PQerrorMessage( this->pDBConn->pPGConn ) ) );
        } // method

        /** @brief Prepare the stmt */
        template <typename... ArgTypes> void prepareStmt( ArgTypes&&... args )
        {
            const size_t uiNumArgs = sizeof...( args ) * this->uiNumRepOfArgs;
            this->iStmtNumParams = (int)uiNumArgs;
            // @fixme getStmtId should be unique across shared libraries...
            this->sStmtName = "STMT_" + std::to_string( (int64_t)&xPG_GLobalEnv ) + "_" +
                              std::to_string( xPG_GLobalEnv.getStmtId( ) );
            // DEBUG: std::cout << "Prepare Stmt " << this->sStmtName << " : " << sStmtText << std::endl;
            vParamValues.resize( uiNumArgs );
            vparamLengths.resize( uiNumArgs );
            vparamFormats.resize( uiNumArgs );

            // Initialize the vInputArgs vector
            for( size_t uiCount = 0; uiCount < uiNumArgs; uiCount++ )
            {
                vInputArgs.emplace_back(
                    StmtArg( vParamValues[ uiCount ], vparamLengths[ uiCount ], vparamFormats[ uiCount ] ) );
            } // for

            argsOIDs( args... );

            // Actually do the prepare call
            if( this->pPGRes ) // avoid memory leak
            {
                PQclear( pPGRes );
                this->pPGRes = NULL;
            } // if

            // Synchronous execution - we wait until completion
            this->pPGRes = PQprepare( this->pDBConn->pPGConn,
                                      sStmtName.c_str( ),
                                      sStmtText.c_str( ),
                                      (int)uiNumArgs,
                                      NULL /* const Oid * paramTypes */ );
            checkResultStatus( );
        } // method

        // for bulk-inserts we require precompiled stmts with the same args repetitively occurring
        size_t uiNumRepOfArgs = 1;

      public:
        /** @brief Constructs a prepared statement using the SQL in rsStmtText
         *  @detail for bulk-inserts we require precompiled stmts with the same args repetitively occurring
         */
        PreparedStmtTmpl( DBPtrType const pDBConn, const std::string& rsStmtText, bool bUseAsyncQuery = false )
            : sStmtName( ), pDBConn( pDBConn ), bUseAsyncQuery( bUseAsyncQuery ), sStmtText( rsStmtText )
        {
            // We must delay the actual creation of the prepared statement until the first parameter binding in order
            // to know the number of parameters.
        } // constructor

        /** @brief Execute the statement after all parameters have been bound successfully.
         *  The existence of parameters has to be initiated via the boolean argument.
         */
        inline size_t exec( bool bMoreThanZeroArgs = true )
        {
            bool bNoArgs = this->iStmtNumParams == 0;
            if( this->pPGRes ) // avoid memory leak
            {
                PQclear( pPGRes );
                this->pPGRes = NULL;
            } // if

            // In async mode we us PQsendQueryPrepared here
            // Asyn Example: https://www.postgresql.org/message-id/20160331195656.17bc0e3b%40slate.meme.com
            if( bUseAsyncQuery )
            {
                // Async operation in single row mode
                if( !PQsendQueryPrepared( this->pDBConn->pPGConn,
                                          this->sStmtName.c_str( ),
                                          this->iStmtNumParams,
                                          bNoArgs ? NULL : &this->vParamValues[ 0 ],
                                          bNoArgs ? NULL : &this->vparamLengths[ 0 ],
                                          bNoArgs ? NULL : &this->vparamFormats[ 0 ],
                                          PG_BINARY_ARG ) ) // results in binary format)
                    throwPostgreSQLError( );
                if( !PQsetSingleRowMode( this->pDBConn->pPGConn ) )
                    throwPostgreSQLError( );
                return 1;
            } // if

            this->pPGRes = PQexecPrepared( this->pDBConn->pPGConn,
                                           this->sStmtName.c_str( ),
                                           this->iStmtNumParams,
                                           bNoArgs ? NULL : &this->vParamValues[ 0 ],
                                           bNoArgs ? NULL : &this->vparamLengths[ 0 ],
                                           bNoArgs ? NULL : &this->vparamFormats[ 0 ],
                                           PG_BINARY_ARG ); // results in binary format
            checkResultStatus( );
            return PQntuples( this->pPGRes );
        } // method

        /** @brief Used in the context of bulk insert stmts for indicating of how many times the args have to occur */
        void setArgsMultiplicator( const size_t uiNumRepeatedOccurences )
        {
            this->uiNumRepOfArgs = uiNumRepeatedOccurences;
        } // method

        /** @brief Binds the arguments args to the parameter block indicated by OFFSET.
         *  Parameter blocks are used for bulk inserts merely. In this case, the SQL-statement arguments are
         *  a repetitive pattern of equal types (the types of a single row). OFFSET N represent represents the
         *  N'th occurrence of such repetitive pattern. The parameters args will be injected (binded) at this N'th
         *  occurrence.
         */
        template <int OFFSET, typename... ArgTypes> inline void bind( ArgTypes&&... args )
        {
            // assert( ( sizeof...( args ) == 0 ) || ( ( iStmtParamCount % (int)( sizeof...( args ) ) == 0 ) &&
            //                                       ( OFFSET * (int)( sizeof...( args ) ) < iStmtParamCount ) ) );
            // stmt name is empty as long as the stmt is not actually prepared
            if( this->sStmtName.empty( ) )
                prepareStmt<ArgTypes&&...>( std::forward<ArgTypes>( args )... );
            bindArgumentsForwarding<OFFSET * sizeof...( args ), ArgTypes&&...>( std::forward<ArgTypes>( args )... );
        } // method

        /** @brief Same as bindStatic, but using a for loop instead of a index sequence.
         *  @details: Using this dynamic approach, it is possible to work with offsets larger than 1000, which bring
         *  most compilers to their limits.
         */
        template <typename... ArgTypes> inline void bindDynamic( int uiOffset, ArgTypes&&... args )
        {
            // assert( ( sizeof...( args ) == 0 ) || ( ( iStmtParamCount % (int)( sizeof...( args ) ) == 0 ) &&
            //                                         ( uiOffset * (int)( sizeof...( args ) ) < iStmtParamCount ) ) );
            // stmt name is empty as long as the stmt is not actually prepared
            if( this->sStmtName.empty( ) )
                prepareStmt<ArgTypes&&...>( std::forward<ArgTypes>( args )... );
            bindArgumentsForwardingRuntime<ArgTypes&&...>( uiOffset * sizeof...( args ),
                                                           std::forward<ArgTypes>( args )... );
        } // method

        /* ArgTypes are the types of the arguments of the query.
         * Return value is the number of rows changed by the statement.
         * Remark: In the current design bind and execute must happen in one method, because the MySQL-bindings
         * happen via references to the method arguments.
         */
        template <typename... ArgTypes> inline uint64_t bindAndExec( ArgTypes&&... args )
        {
            this->bind<0, ArgTypes&&...>( std::forward<ArgTypes>( args )... );
            return this->exec( sizeof...( ArgTypes ) > 0 );
        } // method

        /** @brief Close the statement and free all its resources.
         *  Destructor.
         */
        virtual ~PreparedStmtTmpl( )
        {
            // See https://www.postgresql.org/docs/current/libpq-exec.html for deallocation
            if( this->pPGRes )
            {
                PQclear( pPGRes );
                this->pPGRes = NULL;
            } // if
        } // destructor
    }; // inner class PreparedStmt

    // Standard form of prepared statement, which uses a shared pointer.
    using PreparedStmt = PreparedStmtTmpl<std::shared_ptr<PostgreSQLDBCon>>;

    /** @brief Instance identifies the current owner of a connection the executes a
     *   query in single tuple (row) mode.
     */
    class QuerySingleTupleOwner
    {
      public:
        std::string sOwningStmt;
        QuerySingleTupleOwner( const std::string& rsOwningStmt ) : sOwningStmt( rsOwningStmt )
        {}

        std::string info( )
        {
            return this->sOwningStmt;
        } // method
    }; // inner class

    /** @brief Prepared MySQL query.
     *  @details
     *  Query can be controlled via an optional JSON object:
     *  json{ { MYSQL_EXTRA, { { FLAGS, { NO_RESULT_STORAGE_ON_CLIENT } } } } }
     *  Implementation details can be found at:
     *  https://dev.mysql.com/doc/refman/5.7/en/mysql-stmt-fetch.html
     */
    template <typename DBPtrType, typename... ColTypes> // types of query columns
    class PreparedQueryTmpl : public PreparedStmtTmpl<DBPtrType>
    {
      private:
        static const int NUM_COLS = sizeof...( ColTypes ); // Number of columns of query outcome
        std::tuple<PGRowCell<ColTypes>...> tCellWrappers; // C++ Wrappers for managing the MYSQL
                                                          // binding C-structs (Connection via pointer)
        size_t uiRowCount; // Internal counter that counts the number of rows fetched so far
        int iStatus; // informs about the outcome of the last fetch; initially -1 representing unknown

        /** @brief Performs: 1. Query execution, 2. Binding of results (outcome of query). */
        template <typename... ArgTypes> inline void execBind( ArgTypes... args )
        {
            this->bindAndExec( std::forward<ArgTypes>( args )... ); // served in superclass
            // query result is reachable via this->pPGRes of superclass
            this->bindResult( );
        } // method

        /** @brief Performs: 1. Query execution, 2. binding of results, 3. fetching of first row.
         *  Returns true if a first row exists, false otherwise.
         */
        template <typename... ArgTypes> inline bool execBindFetch( ArgTypes... args )
        {
            this->execBind( std::forward<ArgTypes>( args )... );
            return this->fetchNextRow( );
        } // method

        /** @brief Parses the JSON configuration of he query. */
        void parseJsonConfig( const json& jQueryConfig )
        {
            // Check for the NO_RESULT_STORAGE_ON_CLIENT flag.
#if 0
            if( jQueryConfig.count( MYSQL_EXTRA ) && jQueryConfig[ MYSQL_EXTRA ].count( FLAGS ) )
                for( auto rsString : jQueryConfig[ MYSQL_EXTRA ][ FLAGS ] )
                    if( rsString == NO_RESULT_STORAGE_ON_CLIENT )
                        ; // this->doStoreResult = false;
#endif
        } // method

        /** @brief Translates placeholders from MySQL syntax ('?') to PosgreSQL syntax ('$1', $2' ...)
         */
        std::string adaptPlaceholder( std::string sStmt )
        {
            size_t uiPos = 0;
            size_t uiCounter = 1;
            std::string search = "?";
            while( ( uiPos = sStmt.find( search, uiPos ) ) != std::string::npos )
            {
                std::string sReplace = "$" + std::to_string( uiCounter++ );
                sStmt.replace( uiPos, search.length( ), sReplace );
                uiPos += sReplace.length( );
            } // while
            return sStmt;
        } // method

        friend PostgreSQLDBCon; // give PostgreSQLDBCon access to execBindFetch( )

        const std::shared_ptr<QuerySingleTupleOwner> pxInProgressOwner =
            std::make_shared<QuerySingleTupleOwner>( this->sStmtText );

      public:
        std::tuple<ColTypes...> tCellValues; // tuple that keeps the values of current query row

        /** @brief Constructs a predefined query using the passed SQL statement.
         *  The behavior of the query can be additionally controlled by an optionally passed JSON object.
         *  All placeholders are translated from MySQL syntax ('?') to PosgreSQL syntax ('$1', $2' ...)
         */
        PreparedQueryTmpl( DBPtrType pDBConn, const std::string& rsStmtText, const json& rjConfig = json{} )
            : PreparedStmtTmpl<DBPtrType>(
                  pDBConn, adaptPlaceholder( rsStmtText ), true /* do async */ ), // class superclass constructor
              tCellWrappers( ), // initialized via default constructors (couldn't find better way :-( )
              iStatus( 0 ) // initially 'fetchNextRow not called at all'
        {
            // Connect row cell wrappers and tuple keeping the cell values itself.
            for_each_in_tuple_pairwise(
                tCellWrappers,
                [&]( auto& rFstCell, auto& rSecCell, size_t uiCol ) {
                    // DEBUG: std::cout << "uiColNum:" << uiColNum << " uiCol: " << uiCol << std::endl;
                    rFstCell.init( &rSecCell, uiCol );
                },
                tCellValues );

            // Parse the JSON configuration of the query.
            parseJsonConfig( rjConfig );
        } // constructor

        /** @brief Performs argument binding before fetching. Has to be called after statement execution and
         *  before fetching table rows.
         */
        inline void bindResult( )
        {
            this->uiRowCount = 0; // Reset row counter
            this->iStatus = 0; // status 'fetchNextRow not called at all'
            this->pDBConn->setCommandInProgress(
                pxInProgressOwner ); // Inform DB connection that a query is executed in single tuple mode
        } // method

        /** @brief Fetch next row from server or local client buffer.
         *  Returns true, if a row could be fetched successfully, false otherwise (indicating EOF).
         */
        inline bool fetchNextRow( )
        {
            // former result not used any longer. (avoid memory leak)
            if( this->pPGRes )
            {
                PQclear( this->pPGRes );
                this->pPGRes = NULL;
            } // if

            // If PQgetResult returns null, then there is nothing at all anymore
            if( !( this->pPGRes = PQgetResult( this->pDBConn->pPGConn ) ) )
            {
                this->iStatus = 2; // status 'EOF'
                this->pDBConn->resetCommandInProgress( );
                return false; // done; all results fetched
            } // if

            switch( PQresultStatus( this->pPGRes ) )
            {
                case PGRES_TUPLES_OK: // Indicates next statements results for multiple SQL-stmts
                                      // or end of data for a single stmt.
                    // Don't forget to free each result object with PQclear when done with it.
                    PQclear( this->pPGRes );
                    // PostgreSQL requires that all rows of a result have to be read before executing the next query
                    if( ( this->pPGRes = PQgetResult( this->pDBConn->pPGConn ) ) )
                        throw PostgreSQLError( "PostgreSQL: Logic Error in fetchNextRow. Expected Null pointer." );
                    this->pDBConn->resetCommandInProgress( );
                    this->iStatus = 2; // status 'EOF'
                    return false;

                case PGRES_SINGLE_TUPLE: // Got a row
                    // if( this->uiRowCount == 0 )
                    //     // PQnfields(res) // Check C++ type for query
                    //     ;
                    this->uiRowCount++;

                    // get the actual cell values
                    for_each_in_tuple( tCellWrappers, [&]( auto& rCell ) {
                        if( PQgetisnull( this->pPGRes, 0, (int)rCell.uiColNum ) )
                            rCell.isNull = true;
                        else
                            rCell.store( this->pPGRes /* PQgetvalue( this->pPGRes, 0, (int)rCell.uiColNum ) */ );
                    } ); // for each tuple
                    break;

                default:
                    // Always call PQgetResult until it returns null, even on error
                    this->throwPostgreSQLError( ); // FIXME: Indicate source of error
            } // switch

            this->iStatus = 1; // status 'previous call delivered row'
            return true;
        } // method

        /** @brief Delivers a reference to a tuple holding the data of the current row.
         */
        inline std::tuple<ColTypes...>& getCellValues( )
        {
            if( uiRowCount < 1 )
                throw PostgreSQLError( "PG Logic error. You must first fetch before getting row-data." );
            return this->tCellValues;
        } // method

        /** @brief End Of File. Delivers true if the previous fetch failed. False otherwise.
         *  This kind of eof works in a "look back" fashion; it tells if the previous fetch failed or not.
         */
        inline bool eofLookBack( )
        {
            if( this->iStatus == 0 )
                throw PostgreSQLError(
                    "PG Logic error. eofLookBack() was called without first calling fetchNextRow( )." );
            return this->iStatus == 2;
        } // method
#if 0
        /** @brief End Of File. Delivers true if there are no more rows to fetch.
         *  This kind of eof works in a "look ahead" fashion; it tells if the next fetch will fail or not.
         *  @details Important:Works only if the complete result is stored on client via mysql_stmt_store_result.
         */
        inline bool eofLookAhead( )
        {
            if( !this->doStoreResult )
                throw std::runtime_error( "MySQL: eofLookAhead() requires the call of mysql_stmt_store_result()." );
            return !( uiRowCount < mysql_stmt_num_rows( this->pStmt ) );
        } // method
#endif
        /** @brief explains the execution strategy for the query */
        void explain( )
        {
            // MySQL QueryExplainer<DBPtrType> xQueryExplainer( this->pMySQLDB );
            // MySQL xQueryExplainer.explain( this->sStmtText );
        } // method

        /* Destructor */
        ~PreparedQueryTmpl( )
        {} // destructor
    }; // class ( PreparedQuery )

    template <typename... ColTypes>
    using PreparedQuery = PreparedQueryTmpl<std::shared_ptr<PostgreSQLDBCon>, ColTypes...>;

  public:
    PGconn* pPGConn; // Actual low level PostgreSQL DB Connection
    PGresult* pPGRes; //
    std::string sCurrSchema; //  name of current schema
    // Predefined table existence check
    std::unique_ptr<PreparedQueryTmpl<PostgreSQLDBCon*, int32_t>> pTableExistStmt;

  private:
    /** @brief Check for correct status. In case of difference, throw exception.
     */
    void checkResultStatus( ConnStatusType iStatus )
    {
        if( PQstatus( pPGConn ) != iStatus )
            throw PostgreSQLError( pPGConn );
    } // method

    /** @brief Check for correct result status. In case of difference, throw exception.
     */
    void checkResStatus( ExecStatusType iStatus )
    {
        if( PQresultStatus( pPGRes ) != iStatus )
            throw PostgreSQLError( pPGConn );
    } // method

    /**@brief Make DB connection string
     * FIX ME: For the moment, we always connect to a DB named "postgres"
     */
    std::string makeConnInfo( const json& jConfig )
    {
        // Read configuration from JSON
        std::string sHostName( jConfig.count( HOSTNAME ) > 0 ? jConfig[ HOSTNAME ] : DEFAULT_HOSTNAME );
        std::string sUser( jConfig.count( USER ) > 0 ? jConfig[ USER ] : DEFAULT_USER );
        std::string sPasswd( jConfig.count( PASSWORD ) > 0 ? jConfig[ PASSWORD ] : DEFAULT_PASSWD );
        // unsigned int uiPortNr = jConfig.count( PORT ) > 0 ? jConfig[ PORT ].get<unsigned int>( ) : DEFAULT_PORT;
#if 0
        std::cout << "PostgreSQL::Connection details:\n"
                  << "Hostname: '" << sHostName << "' User: '" << sUser << "' Passwd: '" << sPasswd
                  << "' Port: " << uiPortNr << std::endl;
#endif
        return std::string( "dbname=postgres" ) + " user=" + sUser + " password=" + sPasswd;
    } // method

    /**@brief Opens a MySQL DB connection.
     * @param Connection description in JSON. Default values are used for all missing fields.
     */
    void open( const json& jConfig )
    {
        // Make a connection to the database
        pPGConn = PQconnectdb( makeConnInfo( jConfig ).c_str( ) );
        checkResultStatus( CONNECTION_OK );
    } // method

    std::map<std::string, Oid> mOidMap; // Maps type-names to OID numbers

    /** @brief Read the OID all (name, oid) tuples for the PostgesSQL back-end
     */
    void initOidMap( )
    {
        execSQL( "SELECT typname, oid FROM pg_type", true );
        // currently unsed: int iFields = PQnfields( this->pPGRes );
        for( int i = 0; i < PQntuples( this->pPGRes ); i++ )
            mOidMap.emplace( PQgetvalue( this->pPGRes, i, 0 ), std::stoi( PQgetvalue( this->pPGRes, i, 1 ) ) );
    } // method

    /** @brief Dumps all key-value in the OID map. */
    void dumpOidMap( )
    {
        for( auto& rtKeyValuePair : mOidMap )
            std::cout << rtKeyValuePair.first << "=" << rtKeyValuePair.second << std::endl;
    } // method

#if 0
    /** @brief Sets the directories for CSV file upload during connection startup.
     *  variable tmpdir informs about a temporary directory.
     */
    void setUploadDirs( )
    {
        // TO DO: First check, if there is something in the JSON configuration for the DB.
        // The variable secure_file_priv manages the server side file upload policy:
        //   - NULL           : No server side upload at all.
        //   - "" (empty)	  : Each directory on server side is OK.
        //   - path to folder : The given directory must be used for file upload.
        std::string sSecFilePriv, sTmpDir;
        if( getMySQLVar( "secure_file_priv", sSecFilePriv ) )
        {
            // std::cout << "sSecFilePriv: " << sSecFilePriv << std::endl;
            if( sSecFilePriv != "NULL" )
            {
                if( !sSecFilePriv.empty( ) )
                    pServerDataUploadDir = std::make_unique<fs::path>( sSecFilePriv );
                else
                {
                    if( getMySQLVar( "tmpdir", sTmpDir ) )
                        // ZEUS EXPERIMENT pServerDataUploadDir = std::make_unique<fs::path>( "/mnt/ssd1/tmp" );
                        pServerDataUploadDir = std::make_unique<fs::path>( sTmpDir );
                } // else
            } // if (not NULL)
        } // if (secure_file_priv defined)
    } // method

    // MySQL is missing the "CREATE INDEX IF NOT EXISTS" construct.
    // Therefore, we define a predefined statement for checking the existence of an index.
    std::unique_ptr<PreparedQueryTmpl<MySQLConDB*, int64_t>> pIndexExistStmt;

    std::map<std::string, std::string> mMySQLVars; // local mirror of all MySQL client variables
    unsigned long uiCLientVersion; // version of current client
    std::unique_ptr<fs::path> pServerDataUploadDir; // If null, server upload is not allowed.
    // DEL: std::string sSchemaInUse; // Name of the current schema (managed by useSchema)
#endif

  public:
    PostgreSQLDBCon( const PostgreSQLDBCon& db ) = delete; // no object copies
    PostgreSQLDBCon& operator=( const PostgreSQLDBCon& db ) = delete; // no object assignments

    /** @brief Constructs a MySQL DB connection. Configuration is given via a JSON object */
    PostgreSQLDBCon( const json& jDBConfig = {} ) : pPGConn( NULL ), pPGRes( NULL ), pTableExistStmt( nullptr )
    {
        // See:
        // https://stackoverflow.com/questions/4181951/how-to-check-whether-a-system-is-big-endian-or-little-endian/4181991
        int iVal = 1;
        // big endian endian if true.
        // FIXME: The value storage has to be adapated for big-endian platforms
        if( !( *(char*)&iVal == 1 ) )
            throw std::runtime_error( "The PostgreSQL requires code adaption for your big endian platform  " );

        // Establish connection to database
        open( jDBConfig.count( CONNECTION ) > 0 ? jDBConfig[ CONNECTION ] : json{} );
        // Set always-secure search path, so malicious users can't take control.
        execSQL( "SELECT pg_catalog.set_config('search_path', '', false)", true );
        // Suppress PostgreSQL notices (as e.g. that a schema exists already)
        // https://stackoverflow.com/questions/26150758/suppressing-notice-relation-exists-when-using-create-if-not-exists
        execSQL( "SELECT set_config('client_min_messages', 'error', false)", true );

        initOidMap( );
    } // constructor

    /* Destructor */
    virtual ~PostgreSQLDBCon( )
    {
        do_exception_safe( [&]( ) { this->close( ); } );
    } // destructor

    std::shared_ptr<QuerySingleTupleOwner> pCommInProgrInformer = nullptr;

    /** @brief In PG, fetchNextRow must be called until it delivers false or we get an exception here */
    void setCommandInProgress( std::shared_ptr<QuerySingleTupleOwner> pCommInProgrInformer )
    {
        if( this->pCommInProgrInformer )
            throw PostgreSQLError( "PG Logic Error: Command " + this->pCommInProgrInformer->info( ) + " in progress." );
        this->pCommInProgrInformer = pCommInProgrInformer;
    } // method

    void resetCommandInProgress( )
    {
        this->pCommInProgrInformer = nullptr;
    } // method

    /** @brief: Immediately execute the SQL statement given as argument.
     */
    int execSQL( const std::string& rsSQLStatement, bool bAsQuery = false ) //, bool sSuppressException = false )
    {
        checkDBConn( );
        // Free the memory of the previous exec for avoiding memory leaks.
        if( this->pPGRes )
        {
            PQclear( this->pPGRes );
            this->pPGRes = NULL;
        }

        this->pPGRes = PQexec( pPGConn, rsSQLStatement.c_str( ) );
        checkResStatus( bAsQuery ? PGRES_TUPLES_OK : PGRES_COMMAND_OK );
        return 0;
    } // method

    /** @brief: Closes the database and releases all allocated memory. */
    void close( )
    {
        if( this->pPGRes )
        {
            PQclear( this->pPGRes );
            this->pPGRes = NULL;
        } // if

        if( this->pPGConn )
        {
            // close the connection to the database and cleanup
            PQfinish( this->pPGConn );
            this->pPGConn = NULL;
        } // if
    } // method

    /** @brief Informs the connection about the schema that has to used.
     *  If the schema does not exist already, it is created.
     */
    void useSchema( const std::string& rsSchemaName )
    {
        checkDBConn( );
        execSQL( std::string( "CREATE SCHEMA IF NOT EXISTS " ) + rsSchemaName );
        execSQL( std::string( "SET search_path TO " ) + rsSchemaName + ", public" );
        sCurrSchema = rsSchemaName;
    } // method
#if 0
    /** @brief Delivers the directory that shall be used for data uploads */
    fs::path getDataUploadDir( )
    {
        if( pServerDataUploadDir == nullptr )
            throw MySQLConException( "Server side upload not possible." );
        return *pServerDataUploadDir;
    } // path
#endif
    /** @brief Checks for table existence in current database
     * See: https://stackoverflow.com/questions/20582500/how-to-check-if-a-table-exists-in-a-given-schema
     */
    bool tableExistsInDB( const std::string& sTblName )
    {
        // Precompile the table exist stmt
        if( this->pTableExistStmt == nullptr )
            this->pTableExistStmt = std::make_unique<PreparedQueryTmpl<PostgreSQLDBCon*, int32_t>>(
                this, std::string( "SELECT to_regclass($1)" ) );

        pTableExistStmt->execBindFetch( std::string( this->sCurrSchema ) + "." + sTblName );
        bool bTableExists = !( std::get<0>( this->pTableExistStmt->tCellWrappers ).isNull );
        if( pTableExistStmt->fetchNextRow( ) )
            throw PostgreSQLError( "PostgreSQL: Table exists logic error" );
        return bTableExists;
    } // method

    /* Needs a global lock in concurrent environments */
    bool indexExistsInDB( const std::string& sTblName, const std::string& sIdxName )
    {
        // FIXME if( !( pIndexExistStmt->execBindFetch( sTblName, sIdxName ) ) )
        // FIXME     throw std::runtime_error( "MySQL error: Query in indexExists delivers no result.\n" );
        // FIXME // First element in first row indicates the number of indexes.
        // FIXME return std::get<0>( pIndexExistStmt->tCellValues ) > 0;
        return false;
    } // method

    // FIXME: Change to prepared statements for performance improvements
    inline void startTrxn( )
    {
        // DEBUG: std::cout << "Called startTrxn( ) ... " << std::endl;
        this->execSQL( "BEGIN" );
    } // method

    // FIXME: Change to prepared statements for performance improvements
    inline void commitTrxn( )
    {
        // DEBUG: std::cout << "Called commitTrxn( ) ..." << std::endl;
        this->execSQL( "END" );
    } // method

    /** @brief MySQL does not support WHERE clauses with index creation */
    inline static bool supportsWhereClauseForIndex( )
    {
        return false;
    } // method

    /** @brief Load the file 'sCSVFileName' into the table 'sTblName' */
    void fillTableByFile( const fs::path& sCSVFileName, const std::string sTblName )
    {
        // // Some clients seg. fault if they see the keyword lOCAL
        // bool bUseKeywordLocal = uiCLientVersion != 100141;
        // // Important: Use the generic filename here, because MySQL does not like backslashes in Windows.
        // this->execSQL( std::string( "LOAD DATA " + std::string( bUseKeywordLocal ? "LOCAL " : "" ) + "INFILE \"" +
        //                             sCSVFileName.generic_string( ) + "\" INTO TABLE " + sTblName +
        //                             " FIELDS TERMINATED BY '\t' ENCLOSED BY '\\\'' ESCAPED BY '\\\\'" ) );
    } // method

  protected:
    void checkDBConn( )
    {
        if( pPGConn == NULL )
            throw PostgreSQLError( "PostgreSQL: NULL pointer for DB connection" );
    } // method
}; // class
