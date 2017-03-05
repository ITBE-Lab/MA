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
#include "meta_programming.h"

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

/* Constructs the full filename for a prefix, suffix combination.
 */
std::string fullFileName( const char *pcFileNamePrefix, const char *pcSuffix );

void vWriteStringToFile( const char* pcFileName, std::string &sString );

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

/************************************************************************/
/* Splitter functions for parsing purposes                              */
/************************************************************************/

/* taken from http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
 * There are efficient BOOST solutions for splitting as well.
 * TO DO: use splitAndMap instead of genericSplit
 */
template<class FUNCTION_TYPE>
void genericSplit( const std::string &rsString, const char cDelimiter, FUNCTION_TYPE &&function ) {
	std::string bufferString;
	std::string::const_iterator iterator;

	for ( iterator = rsString.begin(); iterator < rsString.end(); iterator++ ) 
	{
		if ( (const char)*iterator != cDelimiter ) 
		{
			bufferString += *iterator;
		} // if 
		else 
		{
			function( std::move( bufferString) );
			bufferString.clear();
		} // else
	} // for

	if (! bufferString.empty() )
	{
		function( std::move( bufferString) );
	} // if
} // function

/* Taken from http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
 * There are efficient BOOST solutions for splitting as well.
 * Splits the input string according to the given delimiter and applies function to each 
 * text-segment resulting from the splitting.
 * This function is quite inefficient. More efficient would be some modifying (writing '/0') variant with respect to the splited string.
 * However, in this case we have to work with some modifiable buffer.
 * Idea: Implement this by using NonAlteringStringSplittingFunctor.
 */
template<class FUNCTION>
inline void splitAndMap( const std::string &rsSplitString,	// string to split
					     const char cDelimiter,				// delimiter
						 const int iMaxTokenToParse,		// Maximum number of parsed token, deliver 0 for no limitation.
					     FUNCTION &&function				// function applied to each separated text segment.
					   ) 
{
	std::string bufferString; // Use substrings instead of this constructive approach.
	int uiTokenCounter = 0;

	for ( auto iterator = rsSplitString.begin(); iterator < rsSplitString.end();  ) 
	{
		if ( (const char)*iterator != cDelimiter ) 
		{
			bufferString += *iterator;
			iterator++;
		} // if 
		else 
		{
			if ( !( ( iMaxTokenToParse > 0 ) && ( uiTokenCounter + 1 >= iMaxTokenToParse ) ) )
			{
				/* We allow moving away the content of the buffer-string.
				 */
				function( std::move( bufferString ), uiTokenCounter++ );
				
				/* Bring the string back into a known state.
				 * See: http://stackoverflow.com/questions/9168823/reusing-a-moved-container
				 */
				bufferString.clear();
				iterator++;

				/* Move to the first non-delimiter after sequence of delimiter.
				 */
				while( iterator < rsSplitString.end() && (const char)*iterator == cDelimiter ) 
				{
					iterator++;
				} // while
			} // if
			else
			{
				bufferString += *iterator;
				iterator++;
			} // else
		} // else
	} // for

	if (! bufferString.empty() )
	{
		/* We allow moving away the content of the buffer-string.
		 */
		function( std::move( bufferString ), uiTokenCounter );
	} // if
} // function

/* Non-altering string splitter,  i.e. the argument string stays unaltered.
 * Idea: Implement this with the TextSequence call and check, whether we can get significant improvements.
 * bFinalSegmentMayIncludeDelimiter == true : The boolean flag bIsLastField is evaluated and if it inicates true, 
 *											  then we put all remaining text (incl. delimiters) into the field that is marked as last one.
 */
template<bool bFinalSegmentMayIncludeDelimiter>
class NonAlteringStringSplittingFunctor
{
private :
	std::string::const_iterator xStringIterator;	// end of the string to split
	const std::string::const_iterator xStringEnd;	// begin of the split-string
	const char cDelimiter;							// delimiter used for splitting purposes
	std::string sNextToken;							// repeatedly used string buffer
	
	inline void vGetNextToken( bool bIsLastField )
	{
		sNextToken.clear();
		if ( !bFinalSegmentMayIncludeDelimiter || !bIsLastField )
		{
			while ( ( xStringIterator < xStringEnd ) && ( *xStringIterator != cDelimiter ) )
			{
				sNextToken += *(xStringIterator++); // much more efficient than push_back
			} // while
			if ( xStringIterator != xStringEnd )
			{
				xStringIterator++;
			}
		} // if
		else
		{
			while ( xStringIterator < xStringEnd )
			{
				sNextToken += *(xStringIterator++); // much more efficient than push_back
			} // while
		} // else
	} // method

	/* Overloaded getter for an integer value.
	 */
	inline void setAccordingToNextToken( int &rValue, bool bIsLastField )
	{
		/* Extract next substring from split string and call conversion.
		 */
		vGetNextToken( bIsLastField );
		rValue = std::stoi( sNextToken );
	} // method

	/* Overloaded getter for an unsigned long integer value.
	 * Here we catch size_t on 64bit systems. (Hmm... seems to be necessary for Visual C++ 2013)
	 */
	inline void setAccordingToNextToken( unsigned long long &rValue, bool bIsLastField )
	{
		/* Extract next substring from split string and call conversion.
		 */
		vGetNextToken( bIsLastField );
		rValue = std::stoul( sNextToken );
	} // method

	/* Overloaded getter for an unsigned long integer value.
	 */
	inline void setAccordingToNextToken( unsigned long &rValue, bool bIsLastField )
	{
		/* Extract next substring from split string and call conversion.
		 */
		vGetNextToken( bIsLastField );
		rValue = std::stoul( sNextToken );
	} // method

	/* Overloaded getter for a string value.
	 */
	inline void setAccordingToNextToken( std::string &rValue, bool bIsLastField )
	{
		/* Extract next substring from split string and call conversion.
		 */
		vGetNextToken( bIsLastField );
		rValue = std::move( sNextToken );
	} // method

public :
	NonAlteringStringSplittingFunctor
	(	const std::string &rsStringToSplit,
		char cDelimiter
	) : 
		xStringIterator( rsStringToSplit.cbegin() ),
		xStringEnd( rsStringToSplit.cend() ),
		cDelimiter( cDelimiter )
	{} // constructor

	template <typename TupleElementType>
	void operator() ( TupleElementType& rElement, const int I, const int TYPE_SIZE )
	{
		setAccordingToNextToken( rElement, I == ( TYPE_SIZE - 1 ) );
	} // operator

	template <typename TupleElementType>
	void operator() ( TupleElementType& rElement, bool bIsLastElement )
	{
		setAccordingToNextToken( rElement, bIsLastElement );
	} // operator
}; // struct

/* Splits a string into a tuple, where the argument string stays unaltered.
 * The number of types in TupleElementTypes decides the number of applied split operations.
 */
template<typename ...TupleElementTypes>
void vSplitIntoTuple( const std::string &rsStringToSplit,		// string to split by delimiter
					  const char cDelimiter,					// separator used for splitting
					  std::tuple<TupleElementTypes...>& rxTuple	// tuple that receives the outcome of the splitting
				    )
{
	/* NonAlteringStringSplittingFunctor is a functor that overloads the operator() taking 3 arguments.
	 */
	iterateOverTuple( NonAlteringStringSplittingFunctor<true>( rsStringToSplit, cDelimiter ), rxTuple );
}; // function

/* Splits the given string according to the given delimiter and calls the given function with the segments
 * as arguments.
 */
template<typename ...TupleElementTypes, typename Function>
void vSplitAndApply( const std::string &rsStringToSplit,	// string to split by delimiter
					 const char cDelimiter,					// separator used for splitting
					 Function&& function					// function the receives the outcome of the splitting as arguments
				   )
{
	/* First parse into a tuple that works as intermediate storage.
	 */
	std::tuple<TupleElementTypes...> rxTuple;
	iterateOverTuple( NonAlteringStringSplittingFunctor<true>( rsStringToSplit, cDelimiter ), rxTuple );

	/* Call the function with the tuple elements as arguments
	 */
	( TupleUnpackToParameterByReference<TupleElementTypes ...> ( function ) )( rxTuple );
}; // function

/* Splits the given string according to the given delimiter and assigns the values of the separated elements.
 * to args. args are delivered by reference.
 */
template<bool bFinalSegmentMayIncludeDelimiter, typename... Args>
void vSplitAndSet(	const std::string &rsStringToSplit,				// string to split by delimiter
					const char cDelimiter,							// separator used for splitting
					Args&... args									// reference to arguments that shall receive the outcome of the splitting
				 )
{
	metaApplyFunctionToAllArgumentsIndicatingFinalArgument
	(	NonAlteringStringSplittingFunctor<bFinalSegmentMayIncludeDelimiter>( rsStringToSplit, cDelimiter ),	// splitter functor
		std::forward<Args>(args)...																				// arguments iterated over
	);
}; // function

/* Function for range checking.
 * Checks: Whether val is between min and max.
 */
template<typename ParameterType>
void vRangeCheckAndThrowInclusive( const std::string &sText, const ParameterType &xRangeMin, const ParameterType &xVal, const ParameterType &xRangeMax )
{
	if ( xVal < xRangeMin || xVal > xRangeMax )
	{
		throw std::runtime_error(   (((((((std::string( sText ) += "Out of range for value : ") += std::to_string( xVal )) += " range : [ ")
								 += std::to_string( xRangeMin )) += "..") += std::to_string( xRangeMax )) += "]")
								); // runtime error
	} // if
} // template function 

template<typename ParameterType>
void vRangeCheckAndThrowExclusive( const std::string &sText, const ParameterType &xRangeMin, const ParameterType &xVal, const ParameterType &xRangeMax )
{
	if ( xVal < xRangeMin || xVal >= xRangeMax )
	{
		throw std::runtime_error(   (((((((std::string( sText ) += "Out of range for value : ") += std::to_string( xVal )) += " range : [ ")
								 += std::to_string( xRangeMin )) += "..") += std::to_string( xRangeMax )) += ")")
								); // runtime error
	} // if
} // template function 

/* Color related stuff
 */

class RGB_Color
{
public :
	double dRed; // range [0.0 .. 1.0]
	double dGreen; // range [0.0 .. 1.0]
	double dBlue; // range [0.0 .. 1.0]

	/* Delivers the current RGB color as hex-string.
	 */
	std::string sAsHexString()
	{
		uint16_t iRed = std::min( (uint16_t)(dRed * 255), (uint16_t)255 );
		uint16_t iGreen = std::min( (uint16_t)(dGreen * 255), (uint16_t)255 );
		uint16_t iBlue = std::min( (uint16_t)(dBlue * 255), (uint16_t)255 );
		uint32_t iRGBValue = ( (iRed & 0xFF) << 16 ) + ( (iGreen & 0xFF) << 8 ) + (iBlue & 0xFF);
		
		/* Taken from: http://stackoverflow.com/questions/5100718/integer-to-hex-string-in-c
		 */
		std::string sHexString( 6, '0' );
		static const char* digits = "0123456789ABCDEF";
		for ( size_t i = 0, j = ( 6 - 1 ) * 4; i < 6; ++i, j -= 4 )
		{
			sHexString[i] = digits[ ( iRGBValue >> j ) & 0x0F ];
		} // for

		return sHexString;
	} // method

	RGB_Color( double dRed, double dGreen, double  dBlue) :
		dRed( dRed ),
		dGreen( dGreen ),
		dBlue( dBlue )
	{} // constructor
}; // class

RGB_Color GetColour( double dValue, double dMinValue, double dMaxValue );

/* ******************************************
 * Deprecated C I/O functions originating from BWA code
 * ******************************************
 */
#include <zlib.h>

#ifdef __GNUC__
	// Tell GCC to validate printf format string and args
	#define ATTRIBUTE(list) __attribute__ (list)
#else
	#define ATTRIBUTE(list)
#endif

#ifdef _MSC_VER
	#define __func__ __FUNCTION__
#endif

#define xopen(fn, mode) err_xopen_core(__func__, fn, mode)
#define xreopen(fn, mode, fp) err_xreopen_core(__func__, fn, mode, fp)
#define xzopen(fn, mode) err_xzopen_core(__func__, fn, mode)

#define err_fatal_simple(msg) _err_fatal_simple(__func__, msg)
#define err_fatal_simple_core(msg) _err_fatal_simple_core(__func__, msg)

void err_fatal( const char *header, const char *fmt, ... ) ATTRIBUTE( (noreturn) );
void err_fatal_core( const char *header, const char *fmt, ... ) ATTRIBUTE( (noreturn) );
void _err_fatal_simple( const char *func, const char *msg ) ATTRIBUTE( (noreturn) );
void _err_fatal_simple_core( const char *func, const char *msg ) ATTRIBUTE( (noreturn) );
FILE *err_xopen_core( const char *func, const char *fn, const char *mode );
FILE *err_xreopen_core( const char *func, const char *fn, const char *mode, FILE *fp );
gzFile err_xzopen_core( const char *func, const char *fn, const char *mode );
size_t err_fwrite( const void *ptr, size_t size, size_t nmemb, FILE *stream );
size_t err_fread_noeof( void *ptr, size_t size, size_t nmemb, FILE *stream );

int err_gzread( gzFile file, void *ptr, unsigned int len );
int err_fseek( FILE *stream, long offset, int whence );
#define err_rewind(FP) err_fseek((FP), 0, SEEK_SET)
long err_ftell( FILE *stream );
int err_fprintf( FILE *stream, const char *format, ... )
ATTRIBUTE( (format( printf, 2, 3 )) );
int err_printf( const char *format, ... )
ATTRIBUTE( (format( printf, 1, 2 )) );
int err_fputc( int c, FILE *stream );
#define err_putchar(C) err_fputc((C), stdout)
int err_fputs( const char *s, FILE *stream );
#define err_puts(S) err_fputs((S), stdout)
int err_fflush( FILE *stream );
int err_fclose( FILE *stream );
int err_gzclose( gzFile file );

/* Used in the context of the BWA stuff
 */
static inline uint64_t hash_64( uint64_t key )
{
	key += ~(key << 32);
	key ^= (key >> 22);
	key += ~(key << 13);
	key ^= (key >> 8);
	key += (key << 3);
	key ^= (key >> 15);
	key += ~(key << 27);
	key ^= (key >> 31);
	return key;
}