/**
 * @file support.h
 * @brief Implements GzipInputStream, vRangeCheckAndThrowInclusive and vRangeCheckAndThrowExclusive
 * @author Arne Kutzner
 */
#ifdef _MSC_VER
#define NOMINMAX
// #include <windows.h>
#endif

#ifndef SUPPORT_H
#define SUPPORT_H

#include "debug.h"
#include "exported.h"
#include "system.h"

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <sys/stat.h>
#include <vector>
/// @endcond

#ifdef __GNUC__
#include <stdint.h>

// SSE2
#include <emmintrin.h>

#elif _MSC_VER
/* Several data-type definitions of stdint.h not contained in the MSC compiler
 */
typedef unsigned __int8 uint8_t;
typedef signed __int8 int8_t;

typedef unsigned __int32 uint32_t;
typedef signed __int32 int32_t;

typedef unsigned __int16 uint16_t;
typedef signed __int16 int16_t;

typedef unsigned __int64 uint64_t;
#endif

bool EXPORTED fileExists( const std::string& rsFile );

void makeDir( const std::string& rsFile );

/* Constructs the full filename for a prefix, suffix combination.
 */
std::string EXPORTED fullFileName( const char* pcFileNamePrefix, const char* pcSuffix );

/**
 * @brief Function for range checking.
 * @details
 * Checks: Whether val is between min and max.
 */
template <typename ParameterType>
void vRangeCheckAndThrowInclusive( const std::string& sText, const ParameterType& xRangeMin, const ParameterType& xVal,
                                   const ParameterType& xRangeMax )
{
    if( xVal < xRangeMin || xVal > xRangeMax )
    {
        throw std::runtime_error(
            ( ( ( ( ( ( ( std::string( sText ) += "Out of range for value : " ) += std::to_string( xVal ) ) +=
                      " range : [ " ) += std::to_string( xRangeMin ) ) += ".." ) += std::to_string( xRangeMax ) ) +=
              "]" ) ); // runtime error
    } // if
} // template function

/**
 * @brief Function for range checking.
 * @details
 * Checks: Whether val is between min and max.
 */
template <typename ParameterType>
void vRangeCheckAndThrowExclusive( const std::string& sText, const ParameterType& xRangeMin, const ParameterType& xVal,
                                   const ParameterType& xRangeMax )
{
    if( xVal < xRangeMin || xVal >= xRangeMax )
    {
        throw std::runtime_error(
            ( ( ( ( ( ( ( std::string( sText ) += "Out of range for value : " ) += std::to_string( xVal ) ) +=
                      " range : [ " ) += std::to_string( xRangeMin ) ) += ".." ) += std::to_string( xRangeMax ) ) +=
              ")" ) ); // runtime error
    } // if
} // template function

bool EXPORTED is_number( const std::string& s );

/**
 * @brief Loop where the counter value is known during compiletime.
 * @details
 * Example Usage:
 *   template <size_t IDX> struct Exec
 *   {
 *       bool operator( )( int& x )
 *       {
 *           std::cout << x << std::endl;
 *           return true;
 *       } // operator
 *   }; // struct
 *
 *   int main()
 *   {
 *      int x = 10;
 *      bool bComplete = TemplateLoop<3, Exec>::iterate(x);
 *      if(bComplete)
 *          std::cout << "true" << std::endl;
 *      else
 *          std::cout << "false" << std::endl;
 *      return 0;
 *   } // function
 * Prints:
 * 10
 * 10
 * 10
 * true
 *
 *
 * The executed struct can return false in order to break the iteration. If done so iterate returns false.
 * Otherwise iterate returns true.
 */
template <size_t c, template <size_t> class Func> struct TemplateLoop
{
    template <typename... TP_PARAMS> static bool iterate( TP_PARAMS&... rParams )
    {
        if( !TemplateLoop<c - 1, Func>::iterate( rParams... ) )
            return false;
        return Func<c - 1>( )( rParams... );
    } // method
}; // struct

template <template <size_t> class Func> struct TemplateLoop<1, Func>
{
    template <typename... TP_PARAMS> static bool iterate( TP_PARAMS&... rParams )
    {
        return Func<0>( )( rParams... );
    } // method
}; // struct

template <template <size_t> class Func> struct TemplateLoop<0, Func>
{
    template <typename... TP_PARAMS> static bool iterate( TP_PARAMS&... rParams )
    {
        return true;
    } // method
}; // struct

std::string EXPORTED demangle( const char* name );

template <class X> std::string type_name( X* pType, bool bBare = false )
{
    if( pType == nullptr )
        return ( bBare ? "" : "[static type] " ) + demangle( typeid( X ).name( ) );
    return ( bBare ? "" : "[dynamic type] " ) + demangle( typeid( *pType ).name( ) );
} // function

template <class X> std::string type_name( std::shared_ptr<X> pType, bool bBare = false )
{
    return type_name( pType.get( ), bBare );
} // function

template <class X> std::string type_name( bool bBare = false )
{
    return type_name<X>( nullptr, bBare );
} // function

bool EXPORTED ends_with( const std::string& rsX, const std::string& rsEnd );

std::vector<std::string> EXPORTED split( const std::string& sSubject, const std::string sRegex );

/**
 * @brief are we on a big endian system?
 * @details
 * Taken from here: https://stackoverflow.com/questions/1001307/detecting-endianness-programmatically-in-a-c-program.
 * Once we move to C++20 we can use https://en.cppreference.com/w/cpp/types/endian instead.
 */
bool /* EXPORTED */ is_big_endian( );

#endif