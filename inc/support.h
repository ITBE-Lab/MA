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
#include "debug.h"

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

#if 0
/* This can be replaced by std::to_string(...) with C++11
	*/
template <typename TP>
std::string numberToStringDeprecated( TP number )
{
	std::ostringstream ss;
	ss << number;
	return ss.str();
}
/* split 
 * Definition in .cpp file
 * ( Split that alters the input )
 */
void split( std::vector<std::string> &result, std::string str, char delim );

void vDump__m128i( __m128i value );
#endif

/* Constructs the full filename for a prefix, suffix combination.
 */
std::string fullFileName( const char *pcFileNamePrefix, const char *pcSuffix );

#if 0
void vWriteStringToFile( const char* pcFileName, std::string &sString );
#endif
/* The order of the base classes is significant here, because we must initialize filtering_streambuf prior to std::istream
 * According to ISO/IEC 14882:2003(E) section 12.6.2:
 * Then, direct base classes shall be initialized in declaration order as they appear in the base-specifier-list 
 * (regardless of the order of the mem-initializers).
 *
 * If a filtering_streambuf or filtering_stream has mode input, data flows from the chain's end to its beginning.
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

/* An extend form of GzipInputStream for file reading.
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

#if 0

//TODO: this might be worth to keep
/* The counterparts for gzip-file output.
 * Main difference:
 * If it has mode output, data flows in the opposite direction; from the beginning of the chain,
 * through zero or more OutputFilters, towards Sink at the end of the chain.
 */
class GzipOutputStream : boost::iostreams::filtering_streambuf<boost::iostreams::output>,
						 public std::ostream
{
protected :
	void vInitialize( std::ostream &xOutputStream );

public :
	/* The argument of the constructor could be an std::ifstream.
	 * WARNING: The stream must have been opened using the mode std::ios::binary.
	 *          xOutputStream has to exist along with the lifetime of the current object.
	 */
	GzipOutputStream( std::ostream &xOutputStream );

protected :
	/* This constructor is only for derived classes so that these classes can call vInitialize after creating some input stream.
	 */
	GzipOutputStream( );

	virtual ~GzipOutputStream();
}; // class

/* An extend form of GzipOutputStream for file writing.
 * Documentation: See GzipInputFileStream.
 */
class GzipOutputFileStream : public GzipOutputStream
{
private :
	std::ofstream xFileOutputStream;

public :
	GzipOutputFileStream( const char* pcFileName );

	virtual ~GzipOutputFileStream();
}; // class

#endif