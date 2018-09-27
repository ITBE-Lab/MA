/**
 * @file support.h
 * @brief Implements GzipInputStream, vRangeCheckAndThrowInclusive and vRangeCheckAndThrowExclusive
 * @author Arne Kutzner
 */

#ifndef SUPPORT_H
#define SUPPORT_H

#include "util/debug.h"

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <fstream>
#include <iostream>
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

// under gnu EXPORTED is not needed to do anything
#define EXPORTED
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

// under msc we need to export or import a function according to weather we are building the dll or
// using it
#ifdef EXPORT
#define EXPORTED __declspec( dllexport )
#else
#define EXPORTED __declspec( dllimport )
#endif
#endif

bool fileExists( const std::string &rsFile );

void makeDir( const std::string &rsFile );

/* Constructs the full filename for a prefix, suffix combination.
 */
std::string EXPORTED fullFileName( const char *pcFileNamePrefix, const char *pcSuffix );

/**
 * @brief Function for range checking.
 * @details
 * Checks: Whether val is between min and max.
 */
template <typename ParameterType>
void vRangeCheckAndThrowInclusive( const std::string &sText, const ParameterType &xRangeMin,
                                   const ParameterType &xVal, const ParameterType &xRangeMax )
{
    if( xVal < xRangeMin || xVal > xRangeMax )
    {
        throw std::runtime_error( (
            ( ( ( ( ( ( std::string( sText ) += "Out of range for value : " ) += std::to_string(
                          xVal ) ) += " range : [ " ) += std::to_string( xRangeMin ) ) += ".." ) +=
              std::to_string( xRangeMax ) ) += "]" ) ); // runtime error
    } // if
} // template function

/**
 * @brief Function for range checking.
 * @details
 * Checks: Whether val is between min and max.
 */
template <typename ParameterType>
void vRangeCheckAndThrowExclusive( const std::string &sText, const ParameterType &xRangeMin,
                                   const ParameterType &xVal, const ParameterType &xRangeMax )
{
    if( xVal < xRangeMin || xVal >= xRangeMax )
    {
        throw std::runtime_error( (
            ( ( ( ( ( ( std::string( sText ) += "Out of range for value : " ) += std::to_string(
                          xVal ) ) += " range : [ " ) += std::to_string( xRangeMin ) ) += ".." ) +=
              std::to_string( xRangeMax ) ) += ")" ) ); // runtime error
    } // if
} // template function

#endif