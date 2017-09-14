#include <vector>
#include "support.h"

/* Constructs the full file name for some prefix, suffix combination.
 * Returns by value for convenience purposes.
 */
std::string fullFileName( const char *pcFileNamePrefix, const char *pcSuffix )
{
	std::string sFileName( pcFileNamePrefix );
	sFileName.push_back('.');
	return sFileName.append( pcSuffix );
} // method


/* ******************************************
 * GZIP streams on the foundation of boost
 * ******************************************
 */

/* The Implementations for the gzip-input related classes GzipInputStream and GzipInputFileStream
 */
void GzipInputStream::vInitialize( std::istream &xInputStream )
{
	/* We check for the gzip magic.
	 */
	unsigned char aGzipMagic[] = {0, 0};
	xInputStream.read( (char *)aGzipMagic, sizeof( aGzipMagic ) );
	bool bDataInGzipFormat = ( aGzipMagic[0] == 0x1f ) && ( aGzipMagic[1] == 0x8b );

	/* Go back to the beginning of the stream.
	 */
	xInputStream.seekg( 0, xInputStream.beg );
	
	/* If we see the magic at the beginning we insert the appropriate filter into the chain.
	 */
	if ( bDataInGzipFormat )
	{
		push( boost::iostreams::gzip_decompressor() );
	}
	push( xInputStream );
} // private method

GzipInputStream::GzipInputStream( std::istream &xInputStream )
	: boost::iostreams::filtering_streambuf< boost::iostreams::input >(), // 1. initializing filtering_streambuf
	  std::istream( this )				  // 2. Initialize stream and connect it with the boost filter
{
	vInitialize( xInputStream );
} // constructor

GzipInputStream::GzipInputStream( )
	: boost::iostreams::filtering_streambuf< boost::iostreams::input >(), // 1. initializing filtering_streambuf
	  std::istream( this )				  // 2. Initialize stream and connect it with the boost filter
{ } // constructor

GzipInputStream::~GzipInputStream()
{} // virtual destructor

GzipInputFileStream::GzipInputFileStream( const std::string &pcFileName )
	: GzipInputStream( ),	// call the default constructor
	  xFileInputStream( pcFileName, std::ios::in | std::ios::binary )	// open the file input stream
{ 
	vInitialize( xFileInputStream );
} // constructor

bool GzipInputFileStream::is_open()
{
	return xFileInputStream.is_open();
} // method

GzipInputFileStream::~GzipInputFileStream()
{
	xFileInputStream.close();
} // virtual destructor

#if 0
//TODO: this might be worth to keep
/* The Implementations for the gzip-input related classes GzipInputStream and GzipInputFileStream
 */
void GzipOutputStream::vInitialize( std::ostream &xOutputStream )
{
	push( boost::iostreams::gzip_compressor() );
	push( xOutputStream );
} // protected method

GzipOutputStream::GzipOutputStream( std::ostream &xOutputStream )
	: boost::iostreams::filtering_streambuf<boost::iostreams::output>(), // 1. initializing filtering_streambuf
	  std::ostream( this )			 // 2. Initialize stream and connect it with the boost filter
{
	vInitialize( xOutputStream );
} // constructor

GzipOutputStream::GzipOutputStream()
	: boost::iostreams::filtering_streambuf<boost::iostreams::output>(), // 1. initializing filtering_streambuf
	  std::ostream( this )			 // 2. Initialize stream and connect it with the boost filter
	{} // protected default constructor

GzipOutputStream::~GzipOutputStream()
	{} // virtual destructor

GzipOutputFileStream::GzipOutputFileStream( const char* pcFileName )
	: GzipOutputStream( ),
	  xFileOutputStream( pcFileName, std::ios::out | std::ios::binary )
{ 
	vInitialize( xFileOutputStream );
} // constructor

GzipOutputFileStream::~GzipOutputFileStream()
{
	xFileOutputStream.close();
} // virtual destructor

#endif