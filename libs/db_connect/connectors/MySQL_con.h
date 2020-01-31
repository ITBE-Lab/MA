/* Authors: Arne Kutzner and Markus Schmidt
 * Created: Nov. 2019
 * This file is part of the ITBE-Lab code collection.
 * MIT License
 * @file MySQL_con.h
 * @brief Connector for MySQL databases
 */
#pragma once
/* Important Win32/64 Notice:
 * The MySQL headers import winsock2.h.
 * This can create trouble in the context of an additional include of "windows.h".
 * Solutions for this conflict are proposed in:
 * https://stackoverflow.com/questions/5971332/redefinition-errors-in-winsock2-h
 */

// Required for the MySQL headers with MSVC for avoiding conflicts
#ifdef _MSC_VER
#define HAVE_STRUCT_TIMESPEC
#define NOMINMAX
#endif

// Deliver comprehensive error messages
#define VERBOSE_ERROR_MESG

// STL imports
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

// For getting the code working with gcc 6.x.x compiler
#if( __GNUC__ && ( __GNUC__ < 7 ) )
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

#if defined( __GNUC__ ) && ( __GNUC__ < 8 )
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

// Integrate json support
// We must include json here already, because some defines in the MySQL-includes conflict with the json library.
#include <json.hpp>

// MySQL import
#include <mysql.h>
// The header <my_global.h> has been gone since MySQL 8.
// The below typedef keeps the code compiling for older releases of MySQL.
#if !defined( MARIADB_BASE_VERSION ) && !defined( MARIADB_VERSION_ID ) && MYSQL_VERSION_ID >= 80001 &&                 \
    MYSQL_VERSION_ID != 80002
typedef bool my_bool;
#endif

#define DEFAULT_HOSTNAME "localhost"
#define DEFAULT_USER "root"
#define DEFAULT_PASSWD "admin"
#define DEFAULT_DBNAME "defaultDB"
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

// General tips for variadic programming with tuples:
// https://www.murrayc.com/permalink/2015/12/05/modern-c-variadic-template-parameters-and-tuples/

/* Call the function func for each component of the tuple passed as argument.
 * The idea is that func takes its parameter by reference and changes the tuple this way.
 * Taken from: https://stackoverflow.com/questions/26902633/how-to-iterate-over-a-stdtuple-in-c-11/26908596
 * Requires C++14
 * The size is returned in order to suppress an unused-variable warning with GCC
 */
template <class F, class... Ts, std::size_t... Is>
void for_each_in_tuple( std::tuple<Ts...>& tuple, F func, std::index_sequence<Is...> )
{
    auto __attribute__( ( unused ) ) unused = {( func( std::get<Is>( tuple ) ), 0 )...};
} // meta

template <class F, class... Ts> void for_each_in_tuple( std::tuple<Ts...>& tuple, F&& func )
{
    for_each_in_tuple( tuple, std::forward<F>( func ), std::make_index_sequence<sizeof...( Ts )>( ) );
} // meta

/* Variant of for_each_in_tuple, that iterates over two tuples simultaneously.
 * The size is returned in order to suppress an unused-variable warning with GCC
 */
template <class F, class... Ts, std::size_t... Is, class... Tss>
void for_each_in_tuple_pairwise( std::tuple<Ts...>& tuple, F func, std::tuple<Tss...>& tuple1,
                                 std::index_sequence<Is...> )
{
    auto __attribute__( ( unused ) ) unused = {( func( std::get<Is>( tuple ), std::get<Is>( tuple1 ), Is ), 0 )...};
} // meta

template <class F, class... Ts, class... Tss>
void for_each_in_tuple_pairwise( std::tuple<Ts...>& tuple, F&& func, std::tuple<Tss...>& tuple1 )
{
    for_each_in_tuple_pairwise( tuple, std::forward<F>( func ), tuple1, std::make_index_sequence<sizeof...( Ts )>( ) );
} // meta

// Required for workaround with respect to explicit specialization of a template function in a template class.
// See: https://stackoverflow.com/questions/3052579/explicit-specialization-in-non-namespace-scope
template <typename T> struct identity
{
    typedef T type;
}; // struct

/** @brief MySQL Connector Exception
 */
class MySQLConException : public std::runtime_error
{
  private:
    /** @brief Generate a MySQL error message
     */
    std::string makeMessage( MYSQL* pMySQLHandler )
    {
        const char* pErrorText = pMySQLHandler ? mysql_error( pMySQLHandler ) : NULL;

        std::stringstream xBuffer;
        xBuffer << "MySQL database error: " << ( pErrorText ? pErrorText : "UNKNOWN ERROR" ) << std::endl;
        return xBuffer.str( );
    } // method


  public:
    /** @brief MySQLConException exception
     */
    MySQLConException( const std::string& rsErrExplain )
        : std::runtime_error( "MySQL connector error: " + rsErrExplain )
    {}

    /** @brief  Ask MySQL for the error text
     */
    MySQLConException( MYSQL* pMySQLHandler ) : std::runtime_error( makeMessage( pMySQLHandler ) )
    {}
}; // class


/** @brief Executes the functor so that exceptions are swallowed.
 *  @details For destructors so that they do not throw.
 */
template <typename FUNCTOR> inline void do_exception_safe( FUNCTOR&& functor )
{
    try
    {
        functor( );
    } // try
    catch( std::runtime_error& e )
    { // Swallow the exception
        std::cout << e.what( ) << std::endl;
    } // catch
    catch( ... )
    {
        std::cout << "Got unknown exception." << std::endl;
    } // catch
} // function

/** @brief Minimalistic generic Blob (pointer plus size)
 */
struct GenericBlob
{
    char* pBuf; // pointerToBuffer
    size_t uiSize; // size of buffer
}; // struct


/* Basic class for a single cell in row for a query outcome.
 * Full support of C++17 would allow moving all these specializations inside MySQLConDB.
 */
template <typename CellType> class RowCellBase
{

  private:
    /** @brief This class is for cases where we want to have a simple byte buffer without the initialization as it
     * occurs in std::vector.
     */
    class ByteBuffer
    {
      private:
        size_t uiBufSize = 0;
        char* pBuffer = NULL;
        static const size_t SEG_SIZE = 512; // segment size chosen in the context of allocations

      public:
        /** @brief Delivers pointer to buffer. */
        inline char* get( )
        {
            return pBuffer;
        } // method

        /** @brief Returns the size of the actually allocated buffer */
        template <typename Type> inline Type resize( Type uiReq )
        {
            if( uiReq > uiBufSize )
            {
                // uiReq is guaranteed to be one at least
                // allocate multiple of SEG_SIZE
                uiBufSize = ( ( ( uiReq - 1 ) / SEG_SIZE ) + 1 ) * SEG_SIZE;

                // realloc frees automatically in the case of relocation
                pBuffer = static_cast<char*>( realloc( pBuffer, uiBufSize ) );
                if( !pBuffer )
                    throw std::runtime_error( "ByteBuffer reports out of memory for size " +
                                              std::to_string( uiBufSize ) );
            }
            return static_cast<Type>( uiBufSize );
        } // method

        /* Free at destruction */
        ~ByteBuffer( )
        {
            if( pBuffer != NULL )
                free( pBuffer );
        } // destructor
    }; // class

  public:
    MYSQL_BIND* pMySQLBind; // pointer to MySQL bind
    CellType* pCellValue; // pointer to the actual cell value (this pointer refers into tCellValues)
    unsigned long uiLength; // buffer for length in MYSQL_BIND (32 Bit on MSVC, 64 Bit on GCC)
                            // It must be set however for variable-length types, such as BLOBs or STRINGs.
                            // If length is set, mysql_stmt_fetch will write column length into it.
    ByteBuffer pVarLenBuf; // buffer for variable length data (unused for fixed length data)
    my_bool is_null; // indicates if the cell comprises a NULL value
    my_bool error; // 1 if data truncations happened during fetch
    size_t uiColNum; // column number of cell in query
    bool bIsVarSizeCell; // indicates if cell comprises variable sized data

  public:
    /* Initialization of a cell must happen via this init */
    inline void
    init( MYSQL_BIND* pMySQLBind, // pointer to MySQL bind
          CellType* pCellValue, // pointer to actual cell value
          enum_field_types eFieldType, // MYSQL field type of cell
          size_t uiColNum, // column number of cell in query outcome
          bool bIsVarSizeCell = !std::is_arithmetic<CellType>::value ) // MySQL treats arith. data as fixed-size
    {
        memset( pMySQLBind, 0, sizeof( MYSQL_BIND ) ); // clear the MySQL bind record
        this->pMySQLBind = pMySQLBind; // backup pointer to MySQL bind
        this->pCellValue = pCellValue; // backup pointer to cell value
        this->pMySQLBind->buffer_type = eFieldType; // set MySQL datatype indicator
        this->pMySQLBind->length = &this->uiLength; // set location of size buffer
        this->pMySQLBind->is_null = &this->is_null; // set location for null-value indication
        this->pMySQLBind->error = &this->error; // set location for truncation indication
        // Inform MySQL about the cell properties with respect to signed or unsigned.
        // If the MySQL column type and the cell type are different here, MySQL performs truncations.
        this->pMySQLBind->is_unsigned = std::is_unsigned<CellType>::value;

        // The length of the buffer. You don't have to set it for any fixed length buffer: float, double,
        // int, etc. It must be set however for variable-length types, such as BLOBs or STRINGs.
        this->pMySQLBind->buffer_length = 0;

        // For fixed size data we set the MySQL-buffer pointer directly to the cell.
        // Otherwise (for variable size data) we set the MySQL-buffer pointer to the byte buffer.
        this->bIsVarSizeCell = bIsVarSizeCell;
        this->pMySQLBind->buffer = bIsVarSizeCell ? this->pVarLenBuf.get( ) : (void*)pCellValue;

        this->uiColNum = uiColNum; // backup column number
    } // method

    // Dummy method that should never be called.
    inline void storeVarSizeCell( )
    {
        // DEBUG: std::cout << "storeVarSizeCell called for  type: " << typeid( CellType ).name( ) << std::endl;
    } // method
}; // class (RowCellBase)

/* This class is intended to be specialized for specific types. */
template <typename CellType> struct RowCell : public RowCellBase<CellType>
{}; // class

/* The trick here is the overriding of functions in the specialized subclasses.
 * Important: This explicit specializations must occur before any possible implicit specializations or
 *            the compiler complains about "re-specialization".
 */

// boolean value
template <> struct RowCell<bool> : public RowCellBase<bool>
{
    inline void init( MYSQL_BIND* pMySQLBind, bool* pCellValue, size_t uiColNum )
    {
        RowCellBase::init( pMySQLBind, pCellValue, MYSQL_TYPE_TINY, uiColNum );
    } // method
}; // specialized class

// 32 bit integer
template <> struct RowCell<int32_t> : public RowCellBase<int32_t>
{
    inline void init( MYSQL_BIND* pMySQLBind, int32_t* pCellValue, size_t uiColNum )
    {
        RowCellBase::init( pMySQLBind, pCellValue, MYSQL_TYPE_LONG, uiColNum );
    } // method
}; // specialized class

// unsigned 32 bit integer (unsigned info is extracted inside RowCellBase::init )
template <> struct RowCell<uint32_t> : public RowCellBase<uint32_t>
{
    inline void init( MYSQL_BIND* pMySQLBind, uint32_t* pCellValue, size_t uiColNum )
    {
        RowCellBase::init( pMySQLBind, pCellValue, MYSQL_TYPE_LONG, uiColNum );
    } // method
}; // specialized class

// 64 bit integer
template <> struct RowCell<int64_t> : public RowCellBase<int64_t>
{
    inline void init( MYSQL_BIND* pMySQLBind, int64_t* pCellValue, size_t uiColNum )
    {
        RowCellBase::init( pMySQLBind, pCellValue, MYSQL_TYPE_LONGLONG, uiColNum );
    } // method
}; // specialized class

//  unsigned 64 bit integer (unsigned info is extracted inside RowCellBase::init )
template <> struct RowCell<uint64_t> : public RowCellBase<uint64_t>
{
    inline void init( MYSQL_BIND* pMySQLBind, uint64_t* pCellValue, size_t uiColNum )
    {
        RowCellBase::init( pMySQLBind, pCellValue, MYSQL_TYPE_LONGLONG, uiColNum );
    } // method
}; // specialized class

// 64 bit double
template <> struct RowCell<double> : public RowCellBase<double>
{
    inline void init( MYSQL_BIND* pMySQLBind, double* pCellValue, size_t uiColNum )
    {
        RowCellBase::init( pMySQLBind, pCellValue, MYSQL_TYPE_DOUBLE, uiColNum );
    } // method
}; // specialized class

// std::string
template <> struct RowCell<std::string> : public RowCellBase<std::string>
{
    inline void init( MYSQL_BIND* pMySQLBind, std::string* pCellValue, size_t uiColNum )
    {
        RowCellBase::init( pMySQLBind, pCellValue, MYSQL_TYPE_LONG_BLOB, uiColNum );
    } // method

    inline void storeVarSizeCell( )
    {
        pCellValue->assign( this->pVarLenBuf.get( ), this->uiLength );
    } // method
}; // specialized class


/** @brief MqSQL connection to database
 */
class MySQLConDB
{
  public:
    /// @brief Translates C++ to MySQL Column Types
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

      private:
        // list of valid types:
        template <typename Type> static inline std::string getSQLTypeName( identity<Type> );
        // {
        //     return "UNDEFINED_TYPE";
        // } // private method

        static inline std::string getSQLTypeName( identity<std::string> )
        {
            return "LONGTEXT";
        } // private method

        static inline std::string getSQLTypeName( identity<GenericBlob> )
        {
            return "LONGBLOB";
        } // private method

        static inline std::string getSQLTypeName( identity<bool> )
        {
            return "TINYINT";
        } // private method

        static inline std::string getSQLTypeName( identity<int32_t> )
        {
            return "INT";
        } // private method
        static inline std::string getSQLTypeName( identity<uint32_t> )
        {
            return "INT UNSIGNED";
        } // private method

        static inline std::string getSQLTypeName( identity<int64_t> )
        {
            return "BIGINT";
        } // private method

        static inline std::string getSQLTypeName( identity<uint64_t> )
        {
            return "BIGINT UNSIGNED";
        } // private method

        static inline std::string getSQLTypeName( identity<double> )
        {
            return "DOUBLE";
        } // private method
    }; // class
#if 0
    /** @brief Translates C++ to MySQL Column Types
     */
    class TypeTranslator
    {
      public:
        /* Translation of C++ type to SQLite column types.
         * The approach via inlined templates works with older version of GCC as well.
         */
        // uint64_t is too large for sqlite3. Therefore the template
        // template <> std::string SQLTypeName<uint64_t>::name = "INTEGER"
        // this SHALL cause an error for unknown types
        template <typename TP_TYPE> static std::string getSQLTypeName( );

        // list of valid types:
        template <> inline static std::string getSQLTypeName<std::string>( )
        {
            return "LONGTEXT"; // up to 2^32 characters
        } // specialized function
        template <> inline static std::string getSQLTypeName<bool>( )
        {
            return "BOOLEAN";
        } // function

        // numeric:
        // signed
        template <> inline static std::string getSQLTypeName<int8_t>( )
        {
            return "TINYINT";
        } // function
        template <> inline static std::string getSQLTypeName<int16_t>( )
        {
            return "SMALLINT";
        } // function
        template <> inline static std::string getSQLTypeName<int32_t>( )
        {
            return "INT";
        } // function
        template <> inline static std::string getSQLTypeName<int64_t>( )
        {
            return "BIGINT";
        } // function
        // unsigned
        template <> inline static std::string getSQLTypeName<uint8_t>( )
        {
            return "TINYINT";
        } // function
        template <> inline static std::string getSQLTypeName<uint16_t>( )
        {
            return "SMALLINT";
        } // function
        template <> inline static std::string getSQLTypeName<uint32_t>( )
        {
            return "BIGINT";
        } // function

        // floating point:
        template <> inline static std::string getSQLTypeName<double>( )
        {
            return "DOUBLE";
        } // function
        template <> inline static std::string getSQLTypeName<float>( )
        {
            return "FLOAT";
        } // function
    }; // class
#endif
    /** @brief Input argument of a MySQL query
     *  @details
     */
    class StmtArg
    {
      private:
        MYSQL_BIND* const pMySQLBind; // internal MySQL buffer of argument

      public:
        // buffer for length in MYSQL_BIND (32 Bit on MSVC, 64 Bit on GCC)
        // Only required for BLOB and LONGTEXT
        unsigned long uiLength;

        StmtArg( MYSQL_BIND* const pMySQLBind ) : pMySQLBind( pMySQLBind )
        {
            memset( pMySQLBind, 0, sizeof( MYSQL_BIND ) ); // clear the MySQL record
        } // constructor

        /* templated dispatcher.
         * The dispatcher allows using template specializations as replacement for overloading.
         * TIP: This dispatcher can be specialized as well. (See e.g. for NucSeq)
         */
        template <typename Type> inline void set( const Type& dArg )
        {
            // DEBUG: std::cout << "set with dArg type:" << typeid( dArg ).name( ) << std::endl;
            setOverloaded( dArg );
        } // method

      private:
        /** @brief Binder for arguments in the context of statement execution:
         *  @details See: https://dev.mysql.com/doc/refman/5.7/en/c-api-prepared-statement-type-codes.html
         *  Remainder: You have to guarantee the lifetime of the data bound!
         */
        void inline setOverloaded( void* const& p ) = delete; // pointers are not allowed yet

        void inline setOverloaded( const bool& bVal ) // MySQL-type INT
        {
            pMySQLBind->buffer_type = MYSQL_TYPE_TINY; // 8 bit value
            pMySQLBind->buffer = (void*)&bVal;
            pMySQLBind->is_unsigned = false;
        } // method

        void inline setOverloaded( const int32_t& iVal ) // MySQL-type INT
        {
            pMySQLBind->buffer_type = MYSQL_TYPE_LONG; // 32 bit integer
            pMySQLBind->buffer = (void*)&iVal;
            pMySQLBind->is_unsigned = false;
        } // method

        void inline setOverloaded( const uint32_t& iVal ) // MySQL-type INT UNSIGNED
        {
            pMySQLBind->buffer_type = MYSQL_TYPE_LONG; // 32 bit integer
            pMySQLBind->buffer = (void*)&iVal;
            pMySQLBind->is_unsigned = true;
        } // method

        void inline setOverloaded( const int64_t& iVal ) // MySQL-type BIGINT
        {
            pMySQLBind->buffer_type = MYSQL_TYPE_LONGLONG; // 64 bit integer
            pMySQLBind->buffer = (void*)&iVal;
            pMySQLBind->is_unsigned = false;
        } // method

        void inline setOverloaded( const uint64_t& iVal ) // MySQL-type BIGINT
        {
            pMySQLBind->buffer_type = MYSQL_TYPE_LONGLONG; // 64 bit integer
            pMySQLBind->buffer = (void*)&iVal;
            pMySQLBind->is_unsigned = true;
        } // method

        void inline setOverloaded( const double& dVal ) // MySQL-type DOUBLE
        {
            pMySQLBind->buffer_type = MYSQL_TYPE_DOUBLE;
            pMySQLBind->buffer = (void*)&dVal;
        } // method

        // Corresponding head for C-NULL pointer:
        //
        void inline setOverloaded( const std::nullptr_t& p ) // MySQL-type NULL
        {
            pMySQLBind->buffer_type = MYSQL_TYPE_NULL;
            pMySQLBind->buffer = NULL;
        } // method

        // This should work on most architectures.
        // Due to the C++ standard each character of a string is a single byte.
        // Length should deliver the size without the final NULL. (So, we store the string without the final NULL.)
        void inline setOverloaded( const std::string& rsText ) // MySQL-type LONGTEXT
        {
            this->uiLength = static_cast<unsigned long>( rsText.length( ) );
            pMySQLBind->buffer_type = MYSQL_TYPE_LONG_BLOB; // long text
            pMySQLBind->buffer = (void*)rsText.c_str( );
            pMySQLBind->buffer_length = static_cast<unsigned long>( rsText.length( ) );
            pMySQLBind->length = &this->uiLength;
        } // method

        void inline setOverloaded( const GenericBlob& rxBlob ) // MySQL-type LONGBLOB
        {
            this->uiLength = static_cast<unsigned long>( rxBlob.uiSize );
            pMySQLBind->length = &this->uiLength;
            pMySQLBind->buffer_type = MYSQL_TYPE_LONG_BLOB; // long blob
            pMySQLBind->buffer = (void*)rxBlob.pBuf;
            // FIXME: Check for oversized buffers
            pMySQLBind->buffer_length = static_cast<unsigned long>( rxBlob.uiSize );
        } // method
    }; // inner class StmtArg

    /** @brief Prepared MySQL statement
     *  @details Example: https://dev.mysql.com/doc/refman/8.0/en/mysql-stmt-execute.html
     *  DBPtrType should be MySQLConDB*, std::shared_ptr<MySQLConDB> etc.
     */
    template <typename DBPtrType> class PreparedStmtTmpl
    {
      protected:
        /// @brief Generate a detailed MySQL error message specially for statements
        std::string stmtErrMsg( )
        {
            std::string sErrTxt = "MySQL statement error:\n";
#ifdef VERBOSE_ERROR_MESG
            sErrTxt.append( "Affected stmt: " ).append( sStmtText ).append( "\n" );
            sErrTxt.append( "Detailed MySQL error msg:\n" );
#endif
            sErrTxt.append( mysql_stmt_error( pStmt ) ).append( "\n" );
            return sErrTxt;
        } // method

        MYSQL_STMT* pStmt; // pointer to MySQL statement
        std::vector<MYSQL_BIND> vMySQLInpArgBind; // the array of binding records for MySQL itself
        std::vector<StmtArg> vInputArgs; // the array of wrapped binders with references to the MySQL Bind
        DBPtrType const pMySQLDB; // DB Connection that the statement is embedded in
        int iStmtParamCount; // number of parameters of statement
#ifdef VERBOSE_ERROR_MESG
        const std::string sStmtText; // keep the statement text for more comprehensive error messages.
#endif
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

      public:
        /* Constructor */
        PreparedStmtTmpl( DBPtrType const pMySQLDB, const std::string& rsStmtText )
            : pMySQLDB( pMySQLDB ), sStmtText( rsStmtText )
        {
            // DEBUG: std::cout << "Compile statement: " << sStmtText << std::endl;
            pStmt = mysql_stmt_init( this->pMySQLDB->pMySQLHandler );
            if( !pStmt )
                throw std::runtime_error( "mysql_stmt_init(), out of memory" );

            if( mysql_stmt_prepare( pStmt, rsStmtText.c_str( ), static_cast<unsigned long>( rsStmtText.length( ) ) ) )
                throw std::runtime_error( "mysql_stmt_prepare() failed for statement:\n" + rsStmtText +
                                          "\nExplanation:\n" + stmtErrMsg( ) );

            // Get the parameter count from the statement
            iStmtParamCount = mysql_stmt_param_count( pStmt );

            // Adapt the size of the binding vector according the number of parameter
            vMySQLInpArgBind.resize( iStmtParamCount );
            // Create Argument binder for all parameter
            for( int uiCount = 0; uiCount < iStmtParamCount; uiCount++ )
                vInputArgs.push_back( &vMySQLInpArgBind[ uiCount ] );
        } // constructor

        /** @brief Execute the statement after all parameters have been bound successfully.
         */
        inline my_ulonglong exec( void )
        {
            if( mysql_stmt_bind_param( pStmt, &vMySQLInpArgBind[ 0 ] ) )
                throw std::runtime_error( "mysql_stmt_bind_param() failed.\n" + stmtErrMsg( ) );
            // Execute statement
            if( mysql_stmt_execute( pStmt ) )
                throw std::runtime_error( "mysql_stmt_execute, failed\n" + stmtErrMsg( ) );
            // Get the number of affected rows
            // Improvement: Control the actual call via a flag.
            return mysql_stmt_affected_rows( pStmt );
        } // method

        /** @brief Binds the arguments args to the parameter block indicated by OFFSET.
         *  Parameter blocks are used for bulk inserts merely. In this case, the SQL-statement arguments are
         *  a repetitive pattern of equal types (the ypes of a single row). OFFSET N represent represents the
         *  N'th occurrence of such repetitive pattern. The parameter args will be injected (binded) at this N'th
         *  occurrence.
         */
        template <int OFFSET, typename... ArgTypes> inline void bind( ArgTypes&&... args )
        {
            assert( ( sizeof...( args ) == 0 ) || ( ( iStmtParamCount % (int)( sizeof...( args ) ) == 0 ) &&
                                                    ( OFFSET * (int)( sizeof...( args ) ) < iStmtParamCount ) ) );
            bindArgumentsForwarding<OFFSET * sizeof...( args ), ArgTypes&&...>( std::forward<ArgTypes>( args )... );
        } // method

        /* ArgTypes are the types of the arguments of the query.
         * Return value is the number of rows changed by the statement.
         * Remark: In the current design bind and execute must happen in one method, because the MySQL-bindings
         * happen via references to the method arguments.
         */
        template <typename... ArgTypes> inline my_ulonglong bindAndExec( ArgTypes&&... args )
        {
            if( sizeof...( ArgTypes ) != iStmtParamCount )
                throw std::runtime_error( "MySQL - bindAndExec: Mismatch of number of arguments. For statement:\n" +
                                          sStmtText + "\n Actual number: " + std::to_string( sizeof...( ArgTypes ) ) +
                                          " Expected number: " + std::to_string( iStmtParamCount ) );
            this->bind<0, ArgTypes&&...>( std::forward<ArgTypes>( args )... );
            return this->exec( );
            // DEPR // Bind all argument
            // DEPR bindArgumentsForwarding<0, ArgTypes&&...>( std::forward<ArgTypes>( args )... );
            // DEPR if( mysql_stmt_bind_param( pStmt, &vMySQLInpArgBind[ 0 ] ) )
            // DEPR     throw std::runtime_error( "mysql_stmt_bind_param() failed.\n" + stmtErrMsg( ) );
            // DEPR
            // DEPR // Execute statement
            // DEPR if( mysql_stmt_execute( pStmt ) )
            // DEPR     throw std::runtime_error( "mysql_stmt_execute, failed\n" + stmtErrMsg( ) );
            // DEPR
            // DEPR // Get the number of affected rows
            // DEPR // Improvement: Control the actual call via a flag.
            // DEPR return mysql_stmt_affected_rows( pStmt );
        } // method

        /** @brief Close the statement and free all its resources.
         *  Destructor.
         */
        virtual ~PreparedStmtTmpl( )
        {
            if( pStmt )
                if( mysql_stmt_close( pStmt ) )
                    std::cout << "MySQL: Closing of statement failed. Details:" << pMySQLDB->errMsg( ) << std::endl;
        } // destructor
    }; // inner class PreparedStmt

    // Standard form of prepared statement, which uses a shared pointer.
    using PreparedStmt = PreparedStmtTmpl<std::shared_ptr<MySQLConDB>>;

    // HOWTO: https://dev.mysql.com/doc/refman/5.7/en/mysql-stmt-fetch.html
    template <typename DBPtrType, typename... ColTypes> // types of query columns
    class PreparedQueryTmpl : public PreparedStmtTmpl<DBPtrType>
    {
      private:
        static const int NUM_COLS = sizeof...( ColTypes ); // Number of columns of query outcome
        std::array<MYSQL_BIND, NUM_COLS> vMySQLCellBind; // Array of MySQL binding C-structs
                                                         // (represents a row of the outcome of the query)
        std::tuple<RowCell<ColTypes>...> tCellWrappers; // C++ Wrappers for managing the MYSQL
                                                        // binding C-structs (Connection via pointer)
        MYSQL_RES* pPrepareMetaResult; // Pointer to auxiliary C-struct of MySQL
        my_ulonglong uiRowCount; // Internal counter the counts the number of rows fetched so far
        // bool bClientKeepsResult; // true if mysql_stmt_store_result() fetched the result of a query, false otherwise.
        bool doStoreResult; // Use mysql_stmt_store_result() to store the outcome of the query in the client
        int iStatus; // informs about the outcome of the last fetch; initially -1 representing unknown

        /** @brief Performs: 1. Query execution, 2. binding of results.
         */
        template <typename... ArgTypes> inline void execBind( ArgTypes... args )
        {
            this->bindAndExec( std::forward<ArgTypes>( args )... );
            this->bindResult( );
        } // method
        friend MySQLConDB; // give MySQLConDB access to execBindFetch( )

        /** @brief Performs: 1. Query execution, 2. binding of results, 3. fetching of first row.
         *  Returns true if a first row exists, false otherwise.
         */
        template <typename... ArgTypes> inline bool execBindFetch( ArgTypes... args )
        {
            this->execBind( std::forward<ArgTypes>( args )... );
            return this->fetchNextRow( );
        } // method
        friend MySQLConDB; // give MySQLConDB access to execBindFetch( )

      public:
        std::tuple<ColTypes...> tCellValues; // tuple that keeps the values of current query row

        /* Constructor */
        PreparedQueryTmpl( DBPtrType pMySQLDB, const std::string& rsStmtText )
            : PreparedStmtTmpl<DBPtrType>( pMySQLDB,
                                           rsStmtText ), // class superclass constructor
              tCellWrappers( ), // initialized via default constructors (couldn't find better way :-( )
              // bClientKeepsResult( false ) // mysql_stmt_store_result() not called so far
              doStoreResult( true ),
              iStatus( -1 ) // initially unknown
        {
            // Get query related info
            this->pPrepareMetaResult = mysql_stmt_result_metadata( this->pStmt );
            if( !pPrepareMetaResult )
                throw std::runtime_error( "mysql_stmt_result_metadata(), returned no meta information\n" +
                                          this->stmtErrMsg( ) );

            // Get number of columns of query outcome and check correctness
            // IMPROVEMT: Additionally check the individual columns with respect to data-type.
            auto uiColumnCount = mysql_num_fields( this->pPrepareMetaResult );
            if( NUM_COLS != uiColumnCount )
                throw std::runtime_error( "MySQL - The number of columns reported by MySQL does not match the number "
                                          "of template parameters.\nNumber of template para: " +
                                          std::to_string( NUM_COLS ) +
                                          " MySQL column count: " + std::to_string( uiColumnCount ) );

            // Connect the internal MySQL bind-structure with the tuple of row cell wrappers as well as the
            // tuple keeping the cell values itself.
            for_each_in_tuple_pairwise( tCellWrappers,
                                        [&]( auto& rFstCell, auto& rSecCell, size_t uiCol ) {
                                            // std::cout << "uiColNum:" << uiColNum << " uiCol: " << uiCol << std::endl;
                                            rFstCell.init( &vMySQLCellBind[ uiCol ], &rSecCell, uiCol );
                                        },
                                        tCellValues );
        } // constructor

        /** @brief Performs argument binding before fetching. Has to called after statement execution and
         * before fetching table rows.
         */
        inline void bindResult( )
        {
            auto bind = &vMySQLCellBind[ 0 ];
            // Bind the result buffers
            if( mysql_stmt_bind_result( this->pStmt, bind ) )
                throw std::runtime_error( "mysql_stmt_bind_result() failed\n" + this->stmtErrMsg( ) );

            // Fetch the query outcome to the client so that the client delivers all rows without contacting the
            // server anymore. Without this statement each row is fetched separately via server requests.
            // Fetching each row directly from server might be necessary for huge result tables, where
            // the query outcome does not fit into clients memory.
            if( this->doStoreResult )
            {
                if( mysql_stmt_store_result( this->pStmt ) )
                    throw std::runtime_error( "mysql_stmt_store_result() failed\n" + this->stmtErrMsg( ) );
            } // if

            // Reset row counter
            uiRowCount = 0;
        } // method

        /** @brief Fetch next row from server or local client buffer.
         *  Returns true, if a row could be fetched successfully, false otherwise (indicating EOF).
         */
        inline bool fetchNextRow( )
        {
            // Fetch next row
            this->iStatus = mysql_stmt_fetch( this->pStmt );

            if( this->iStatus == 1 )
                // Something went seriously wrong.
                throw std::runtime_error( "mysql_stmt_fetch() failed.\n" + this->stmtErrMsg( ) );

#ifdef REPORT_ROWS_WITH_NULL_VALUES
            // Check for possible NULL values
            // Currently, there is no solution for these cases.
            for_each_in_tuple( tCellWrappers, [&]( auto& rCell ) {
                if( rCell.is_null )
                    throw std::runtime_error( "fetchNextRow reports a cell with NULL value." );
            } ); // for each tuple
#endif
            // Truncations occur in the case of undersized buffers.
            if( this->iStatus == MYSQL_DATA_TRUNCATED )
            {
                // Find the column(s) responsible and resize the buffer appropriately.
                for_each_in_tuple( tCellWrappers, [&]( auto& rCell ) {
                    if( rCell.error ) // error identifies the cell(s), where we have truncation
                    {
                        // Truncation can occur for variable size data only.
                        if( !rCell.bIsVarSizeCell )
                            throw std::runtime_error( "fetchNextRow reports truncation for fixed size data." );

                        // Resize the cell buffer appropriately and update the buffer pointer
                        rCell.pMySQLBind->buffer_length = rCell.pVarLenBuf.resize( rCell.uiLength );
                        rCell.pMySQLBind->buffer = rCell.pVarLenBuf.get( );

                        // Fetch the truncated data in complete form
                        // Possible optimization: Use additionally the offset parameter for avoiding repeated fetching
                        if( mysql_stmt_fetch_column( this->pStmt, rCell.pMySQLBind,
                                                     static_cast<unsigned int>( rCell.uiColNum ), 0 ) )
                            throw std::runtime_error( "mysql_stmt_fetch_column() failed\n" + this->stmtErrMsg( ) );
                    } // if
                } ); // for each cell in tuple

                // Inform MySQL with respect to buffer size changes by updating binding infos.
                if( mysql_stmt_bind_result( this->pStmt, &vMySQLCellBind[ 0 ] ) )
                    throw std::runtime_error( "mysql_stmt_bind_result() failed\n" + this->stmtErrMsg( ) );

                this->iStatus = 0; // everything fine now ...
            } // if

            // Pick variable size data (Text or Blob) from cell buffers.
            if( this->iStatus == 0 )
            {
                for_each_in_tuple( tCellWrappers, [&]( auto& rCell ) {
                    if( rCell.bIsVarSizeCell )
                        rCell.storeVarSizeCell( );
                } ); // for each tuple

                uiRowCount++; // increment counter for number of rows fetched
            } // if

            return this->iStatus == 0; // iStatus is 0 (OK) or MYSQL_NO_DATA now
        } // method

        /** @brief Delivers a reference to a tuple holding the data of the current row.
         */
        inline std::tuple<ColTypes...>& getCellValues( )
        {
            if( uiRowCount < 1 )
                throw std::runtime_error( "MySQL: Logic error. You must first fetch before getting row-data." );
            if( this->iStatus != 0 )
                throw std::runtime_error(
                    "MySQL: Logic error. You tried to receive a row although the previous fetch failed." );
            return this->tCellValues;
        } // method

        /** @brief End Of File. Delivers true if the previous fetch failed. False otherwise.
         *  This kind of eof works in a "look back" fashion; it tells if the previous fetch failed or not.
         */
        inline bool eofLookBack( )
        {
            if( this->iStatus == -1 )
                throw std::runtime_error(
                    "MySQL: Logic error. eofLookBack() was called without first calling fetchNextRow( )." );
            return this->iStatus != 0;
        } // method

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

        /* Destructor */
        ~PreparedQueryTmpl( )
        {
            // Free the prepared result metadata
            if( this->pPrepareMetaResult )
                mysql_free_result( pPrepareMetaResult );
            // Release memory associated with the result set produced by execution of the prepared statement.
            // if( bClientKeepsResult )
            //     if( mysql_stmt_free_result( this->pStmt ) )
            //         // FIXME: Write the error to some form of log file
            //         std::cout << ( "mysql_stmt_free_result() failed in destructor.\n" + this->stmtErrMsg( ) )
            //                   << std::endl;
        } // destructor
    }; // inner class (PreparedQuery)

    template <typename... ColTypes> using PreparedQuery = PreparedQueryTmpl<std::shared_ptr<MySQLConDB>, ColTypes...>;

  private:
    MYSQL* pMySQLHandler; // MySQL handler belonging to the current connection.

    /** @brief Generate a detailed MySQL error message
     */
    std::string errMsg( )
    {
        const char* pErrorText = pMySQLHandler ? mysql_error( pMySQLHandler ) : NULL;

        std::stringstream xBuffer;
        xBuffer << "MySQL error: " << ( pErrorText ? pErrorText : "UNKNOWN ERROR" ) << std::endl;
        return xBuffer.str( );
    } // method

    /** @param rsHostName - can be either a host name or an IP address. Passing the NULL value or the string
     * "localhost" to this parameter, the local host is assumed. When possible, pipes will be used instead of the
     * TCP/IP protocol.
     */
    void open( const std::string& rsHostName, const std::string& rsUser, const std::string& rsPasswd,
               const std::string& rsDBName, unsigned int uiPortNr )
    {
        checkDBCon( );
        if( !( mysql_real_connect( pMySQLHandler, rsHostName.c_str( ), rsUser.c_str( ), rsPasswd.c_str( ), NULL,
                                   uiPortNr, NULL, 0 ) ) )
            throw MySQLConException( pMySQLHandler );
    } // method

    /** @brief Tells the connection about the schema that has to used. */
    void setupDB( const std::string& rsDBName )
    {
        checkDBCon( );
        if( mysql_select_db( pMySQLHandler, rsDBName.c_str( ) ) )
            // Create the database if not existing.
            execSQL( "CREATE DATABASE " + rsDBName );

        // Tell MySQL to USE this database
        execSQL( "USE " + rsDBName );
    } // method

    /** @brief Mirrors all MySQL variables locally be executing the "SHOW VARIABLES" statement and storing the outcome
     *  in the std::map mMySQLVars.
     */
    void mirrorMySQLVars( )
    {
        auto pShowVariablesStmt =
            std::make_unique<PreparedQueryTmpl<MySQLConDB*, std::string, std::string>>( this, "SHOW VARIABLES" );
        pShowVariablesStmt->execBind( );
        while( pShowVariablesStmt->fetchNextRow( ) )
            mMySQLVars.emplace( std::get<0>( pShowVariablesStmt->getCellValues( ) ),
                                std::get<1>( pShowVariablesStmt->getCellValues( ) ) );
    } // method

    /** @brief Casts a MySQL (in rsMySQLVal) string to a C++ string (in rCPPVal) */
    void MySQLVarCast( const std::string& rsMySQLVal, std::string& rCPPVal )
    {
        rCPPVal = rsMySQLVal;
    } // method

    /** @brief Dumps all key-value for mirrored MYSQL vars. */
    void dumpMirroredMySQLVars( )
    {
        for( auto& rtKeyValuePair : mMySQLVars )
            std::cout << rtKeyValuePair.first << "=" << rtKeyValuePair.second << std::endl;
    } // method

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
            std::cout << "sSecFilePriv: " << sSecFilePriv << std::endl;
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

  public:
    MySQLConDB( const MySQLConDB& db ) = delete; // no object copies
    MySQLConDB& operator=( const MySQLConDB& db ) = delete; // no object assignments

    /* Parameterized constructor */
    MySQLConDB( const std::string& rsHostName, const std::string& rsUser, const std::string& rsPasswd,
                const std::string& rsDBName, unsigned int uiPortNr )
        : pMySQLHandler( NULL ), pServerDataUploadDir( nullptr )
    {
        // Initialize the connector
        if( !( pMySQLHandler = mysql_init( NULL ) ) )
            throw MySQLConException( "Initialization of MySQL connection failed" );

        // Establish connection to database
        this->open( rsHostName, rsUser, rsPasswd, rsDBName, uiPortNr );

        // Use the database with name rsDBName. If it does not exists, create it.
        this->setupDB( rsDBName );

        // Compile the statement that checks for the existence of an index.
        // See: https://dba.stackexchange.com/questions/24531/mysql-create-index-if-not-exists
        this->pIndexExistStmt = std::make_unique<PreparedQueryTmpl<MySQLConDB*, int64_t>>(
            this, "SELECT COUNT( 1 ) IndexIsThere FROM INFORMATION_SCHEMA.STATISTICS "
                  "WHERE table_schema = DATABASE() AND table_name = ? AND index_name = ?" );

        // Mirror all client side MySQL variables locally.
        this->mirrorMySQLVars( );

        // Collect client info in numeric form.
        this->uiCLientVersion = mysql_get_client_version( );
        std::cout << "MySQL client version: " << uiCLientVersion << std::endl;

        // Find out the appropriate directory for CSV data uploads.
        this->setUploadDirs( );
    } // constructor

    // /* Default constructor */
    // MySQLConDB( ) : MySQLConDB( DEFAULT_HOSTNAME, DEFAULT_USER, DEFAULT_PASSWD, DEFAULT_DBNAME, DEFAULT_PORT )
    // {} // constructor

    /* Constructor that allows specifying a database name. */
    MySQLConDB( const std::string& rsDBName = "" )
        : MySQLConDB( DEFAULT_HOSTNAME, DEFAULT_USER, DEFAULT_PASSWD, rsDBName.empty( ) ? DEFAULT_DBNAME : rsDBName,
                      DEFAULT_PORT )
    {} // constructor

    /* Destructor */
    virtual ~MySQLConDB( )
    {
        do_exception_safe( [&]( ) { this->close( ); } );

        // For really freeing all memory allocated by MySQL ...
        // See: https://stackoverflow.com/questions/8554585/mysql-c-api-memory-leak
        mysql_library_end( );
    } // destructor

    /** @brief: Immediately execute the SQL statement giver as argument.
     *  Must be followed by mysql_store_result or mysql_use_result() for receiving data from DB.
     */
    int execSQL( const std::string& rsSQLStatement, bool sSuppressException = false )
    {
        checkDBCon( );
        auto iErrorCode = mysql_query( pMySQLHandler, rsSQLStatement.c_str( ) );
        if( iErrorCode && !sSuppressException )
            throw MySQLConException( pMySQLHandler );
        return iErrorCode;
    } // method

    /** @brief: Closes the database and releases all allocated memory. */
    void close( )
    {
        if( pMySQLHandler )
        {
            // releases all allocated resources
            mysql_close( pMySQLHandler );
            pMySQLHandler = NULL;
        } // if
    } // method
      // Quick hack
    // https://stackoverflow.com/questions/32737478/how-should-i-tackle-secure-file-priv-in-mysql
    // AND: https://forums.mysql.com/read.php?152,674208,674208
    /** @brief Reads a MySQL variable and returns it value in rVal.
     *  Returns true if the reading succeeded, false otherwise.
     */
    template <typename TypeOfVar> bool getMySQLVar( const std::string rsVarName, TypeOfVar& rVal )
    {
        auto pSearch = this->mMySQLVars.find( rsVarName );
        if( pSearch != this->mMySQLVars.end( ) )
        {
            MySQLVarCast( pSearch->second, rVal );
            return true; // search succeeded
        } // if
        else
            return false; // search failed
    } // method

    fs::path getDataUploadDir( )
    {
        if( pServerDataUploadDir == nullptr )
            throw MySQLConException( "Server side upload not possible." );
        return *pServerDataUploadDir;
    } // path

    /** @brief Checks for table existence in current database */
    bool tableExistsInDB( const std::string& sTableName )
    {
        auto iErrCode = execSQL( "SELECT 1 FROM " + sTableName + " LIMIT 1", // triggers error, if table does not exist
                                 true ); // suppress exception
        if( !iErrCode )
        {
            // We must get the result for avoiding later sync-errors
            MYSQL_RES* __attribute__( ( unused ) ) result = mysql_store_result( pMySQLHandler );
        } // if
        return !iErrCode;
    } // method

    /* Needs a global lock in concurrent environments */
    bool indexExistsInDB( const std::string& sTblName, const std::string& sIdxName )
    {
        if( !( pIndexExistStmt->execBindFetch( sTblName, sIdxName ) ) )
            throw std::runtime_error( "MySQL error: Query in indexExists delivers no result.\n" );
        // First element in first row indicates the number of indexes.
        return std::get<0>( pIndexExistStmt->tCellValues ) > 0;
    } // method

    /** @brief Delivers the primary key of the last inserted row. (On DB level, not table level!) */
    my_ulonglong getLastAutoIncrVal( ) // my_ulonglong == uint64_t
    {
        // See: https://dev.mysql.com/doc/refman/8.0/en/mysql-insert-id.html
        my_ulonglong uiLastAutoIncrVal = mysql_insert_id( pMySQLHandler );
        if( !uiLastAutoIncrVal )
            throw std::runtime_error( "MySQL error:getLastAutoIncrVal() failed (returned 0)\n" );
        return mysql_insert_id( pMySQLHandler );
    } // method

    // FIXME: Change to prepared statements for performance improvements
    inline void startTrxn( )
    {
        // DEBUB: std::cout << "Called startTrxn( ) ... " << std::endl;
        this->execSQL( "SET autocommit=0" );
        this->execSQL( "SET unique_checks=0" );
        this->execSQL( "SET foreign_key_checks=0" );
    } // method

    // FIXME: Change to prepared statements for performance improvements
    inline void commitTrxn( )
    {
        // DEBUG: std::cout << "Called commitTrxn( ) ..." << std::endl;
        this->execSQL( "COMMIT" );
        this->execSQL( "SET unique_checks=1" );
        this->execSQL( "SET foreign_key_checks=1" );
    } // method

    /** @brief MySQL does not support WHERE clauses with index creation */
    inline static bool supportsWhereClauseForIndex( )
    {
        return false;
    } // method

    /** @brief Load the file 'sCSVFileName' into the table 'sTblName' */
    void fillTableByFile( const fs::path& sCSVFileName, const std::string sTblName )
    {
        // Some clients seg. fault if they see local
        bool bUseKeywordLocal = uiCLientVersion != 100141;
        // Important: Use the generic filename here, because MySQL does not like backslashes in Windows.
        this->execSQL( std::string( "LOAD DATA " + std::string( bUseKeywordLocal ? "LOCAL " : "" ) + "INFILE \"" +
                                    sCSVFileName.generic_string( ) + "\" INTO TABLE " + sTblName ) );
    } // method

  protected:
    void checkDBCon( )
    {
        if( pMySQLHandler == NULL )
            throw MySQLConException( "No valid handler reference" );
    } // method
}; // class
