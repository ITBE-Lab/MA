/** 
 * @file support.h
 * @brief Implements GzipInputStream, vRangeCheckAndThrowInclusive and vRangeCheckAndThrowExclusive
 * @author Arne Kutzner
 */
#pragma once

#include <string>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#include <fstream>
#include "util/debug.h"

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


/* Constructs the full filename for a prefix, suffix combination.
 */
std::string fullFileName( const char *pcFileNamePrefix, const char *pcSuffix );

/**
 * @brief reads gzip files from disk
 * @details
 * The order of the base classes is significant here,
 * because we must initialize filtering_streambuf prior to std::istream
 * According to ISO/IEC 14882:2003(E) section 12.6.2:
 * Then, direct base classes shall be initialized in 
 * declaration order as they appear in the base-specifier-list 
 * (regardless of the order of the mem-initializers).
 *
 * If a filtering_streambuf or filtering_stream has mode input,
 * data flows from the chain's end to its beginning.
 * So, here the data flow towards the current (objects) istream.
 */
class GzipInputStream : protected boost::iostreams::filtering_streambuf<boost::iostreams::input>,
                        public std::istream
{
protected :
    /* Initializes the gzip filter. Look for the gzip-magic if this absent we work in some uncompressed mode.
     * In derived classes this method together with the default constructor can be used for delayed stream connecting.
     */
    void vInitialize( std::istream &xInputStream );

public :
    /* The argument of the constructor could be an std::ifstream.
     * WARNING: The stream must have been opened using the mode std::ios::binary.
     *          xInputStream has to exist along with the lifetime of the current object.
     */
    GzipInputStream( std::istream &xInputStream );

protected :
    /* This constructor is only for derived classes so that these classes can call vInitialize after creating some input stream.
     */
    GzipInputStream( );

public :
    virtual ~GzipInputStream();
}; // class

/**
 * @brief An extend form of GzipInputStream for file reading.
 */
class GzipInputFileStream : public GzipInputStream
{
private :
    std::ifstream xFileInputStream;

public :
    GzipInputFileStream( const std::string &pcFileName );

    bool is_open();

    virtual ~GzipInputFileStream();
}; // class

/**
 * @brief Function for range checking.
 * @details
 * Checks: Whether val is between min and max.
 */
template<typename ParameterType>
void vRangeCheckAndThrowInclusive(
        const std::string &sText,
        const ParameterType &xRangeMin,
        const ParameterType &xVal,
        const ParameterType &xRangeMax
    )
{
       if ( xVal < xRangeMin || xVal > xRangeMax )
       {
               throw std::runtime_error(   (((((((std::string( sText ) += "Out of range for value : ") += std::to_string( xVal )) += " range : [ ") += std::to_string( xRangeMin )) += "..") += std::to_string( xRangeMax )) += "]") ); // runtime error
       } // if
} // template function

/**
 * @brief Function for range checking.
 * @details
 * Checks: Whether val is between min and max.
 */
template<typename ParameterType>
void vRangeCheckAndThrowExclusive( const std::string &sText, const ParameterType &xRangeMin, const ParameterType &xVal, const ParameterType &xRangeMax )
{
       if ( xVal < xRangeMin || xVal >= xRangeMax )
       {
               throw std::runtime_error(   (((((((std::string( sText ) += "Out of range for value : ") += std::to_string( xVal )) += " range : [ ") += std::to_string( xRangeMin )) += "..") += std::to_string( xRangeMax )) += ")") ); // runtime error
       } // if
} // template function