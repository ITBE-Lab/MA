/** 
 * @file support.cpp
 * @author Arne Kutzner
 */
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>

#include "util/support.h"
#include "util/exception.h"
/* Constructs the full file name for some prefix, suffix combination.
 * Returns by value for convenience purposes.
 */
std::string fullFileName( const char *pcFileNamePrefix, const char *pcSuffix )
{
    std::string sFileName( pcFileNamePrefix );
    sFileName.push_back('.');
    return sFileName.append( pcSuffix );
} // method


#if USE_BOOST_GZIP == ( 1 )
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
    
#if USE_BOOST_GZIP == ( 1 )
    /* If we see the magic at the beginning we insert the appropriate filter into the chain.
     */
    if ( bDataInGzipFormat )
    {
        push( boost::iostreams::gzip_decompressor() );
    }// if
    push( xInputStream );
#else
    if ( bDataInGzipFormat )
    {
        throw fasta_reader_exception(
                "Trying to read Gzipped file; but did not compile with boost gzip enabled"
            );
    }// if
#endif
} // private method

GzipInputStream::GzipInputStream( std::istream &xInputStream )
    : boost::iostreams::filtering_streambuf< boost::iostreams::input >(), // 1. initializing filtering_streambuf
      std::istream( this )                  // 2. Initialize stream and connect it with the boost filter
{
    vInitialize( xInputStream );
} // constructor

GzipInputStream::GzipInputStream( )
    : boost::iostreams::filtering_streambuf< boost::iostreams::input >(), // 1. initializing filtering_streambuf
      std::istream( this )                  // 2. Initialize stream and connect it with the boost filter
{ } // constructor

GzipInputStream::~GzipInputStream()
{} // virtual destructor

GzipInputFileStream::GzipInputFileStream( const std::string &pcFileName )
    : GzipInputStream( ),    // call the default constructor
      xFileInputStream( pcFileName, std::ios::in | std::ios::binary )    // open the file input stream
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
#endif

bool fileExists(const std::string& rsFile)
{
    struct stat buffer;
    return (stat (rsFile.c_str(), &buffer) == 0);
}// function

void makeDir(const std::string& rsFile)
{
    mode_t nMode = 0733; // UNIX style permissions
    int nError = 0;
    #if defined(_WIN32)
        nError = _mkdir(rsFile.c_str()); // can be used on Windows
    #else 
        nError = mkdir(rsFile.c_str(), nMode); // can be used on non-Windows
    #endif
    if (nError != 0) {
        throw NullPointerException(std::string("Could not create Dir: ").append(rsFile).c_str());
    }// if
}// function
