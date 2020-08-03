/* Authors: Arne Kutzner and Markus Schmidt
 * Created: July 2020
 * This file is part of the ITBE-Lab code collection.
 * MIT License
 * @db_common.h
 * @brief General supportive stuff in the context of DB support
 */

 // Pretty printing for tables
#include <assert.h>
#include <fort.hpp>

 /* If we're not using GNU C, elide __attribute__ */
#ifndef __GNUC__
#define __attribute__( x ) /* NOTHING */
#endif 

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
    auto __attribute__( ( unused ) ) unused = { ( func( std::get<Is>( tuple ) ), 0 )... };
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
    auto __attribute__( ( unused ) ) unused = { ( func( std::get<Is>( tuple ), std::get<Is>( tuple1 ), Is ), 0 )... };
} // meta

template <class F, class... Ts, class... Tss>
void for_each_in_tuple_pairwise( std::tuple<Ts...>& tuple, F&& func, std::tuple<Tss...>& tuple1 )
{
    for_each_in_tuple_pairwise( tuple, std::forward<F>( func ), tuple1, std::make_index_sequence<sizeof...( Ts )>( ) );
} // meta


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


/** @brief This class is for cases where we want to have a simple byte buffer without the initialization as it
 *  occurs in std::vector.
 */
class ByteBuffer
{
  private:
    size_t uiBufSize = 0;
    char* pBuffer = NULL;
    static const size_t SEG_SIZE = 512; // segment size chosen in the context of allocations

  public:
    /** @brief Delivers pointer to buffer. */
    inline char* get( ) const
    {
        return pBuffer;
    } // method

    /** @brief Returns the size of the actually allocated buffer */
    template <typename Type> inline Type resize( Type uiReq )
    {
        if( uiReq > uiBufSize )
        {
            // uiReq is guaranteed to be greater equal one
            // allocate multiple of SEG_SIZE
            uiBufSize = ( ( ( uiReq - 1 ) / SEG_SIZE ) + 1 ) * SEG_SIZE;

            // realloc frees automatically in the case of relocation
            pBuffer = static_cast<char*>( realloc( pBuffer, uiBufSize ) );
            if( !pBuffer )
                throw std::runtime_error( "ByteBuffer reports out of memory for size " + std::to_string( uiBufSize ) );
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


/** @brief Simple table, where all cells keep text. */
class SimpleTextTable
{
    std::vector<std::string> vColHeads;
    std::vector<std::vector<std::string>> vTblBody;

  public:
    /** @brief Constructs a simple text table. */
    SimpleTextTable( const std::vector<std::string>& rvColHeads ) : vColHeads( rvColHeads )
    {}

    /** @brief Add a fresh row to the table. */
    void addRow( const std::vector<std::string>& rvRow )
    {
        assert( rvRow.size( ) == vColHeads.size( ) );
        vTblBody.push_back( rvRow );
    } // method

    /** @brief Print the table in formated layout. */
    void print( )
    {
        fort::char_table xTable;
        xTable << fort::header;
        // print table head
        for( auto& rsStr : vColHeads )
            xTable << rsStr;
        xTable << fort::endr;
        // print table rows
        for( auto& rvRow : vTblBody )
        {
            for( auto& rsStr : rvRow )
                xTable << rsStr;
            xTable << fort::endr;
        } // for
        std::cout << xTable.to_string( ) << std::endl;
    } // method
}; // class ( SimpleTextTable )
