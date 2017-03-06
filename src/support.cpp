#include <vector>
#include "support.h"


/* Splits rsString according to cDelimiter and stores th outcome in rResultVector.
 */
void split( std::vector<std::string> &rResultVector, const std::string &rsString, const char cDelimiter )
{
	rResultVector.clear();
	genericSplit
	(	rsString, cDelimiter,
		[&rResultVector]( std::string &&rsBufferString )
		{
			rResultVector.emplace_back( rsBufferString );
		} // lambda
	); // function call
} // function

#if 0 // deprecated
/* taken from http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
 * There are efficient BOOST solutions for splitting as well.
 */
void split( std::vector<std::string> &rResultVector, const std::string &rsString, const char cDelimiter ) {
	std::string bufferString;
	std::string::const_iterator iterator;
	rResultVector.clear();

	for ( iterator = rsString.begin(); iterator < rsString.end(); iterator++ ) 
	{
		if ( (const char)*iterator != cDelimiter ) 
		{
			bufferString += *iterator;
		} // if 
		else 
		{
			rResultVector.emplace_back(bufferString);
			bufferString.clear();
		} // else
	} // for

	if (! bufferString.empty() )
	{
		rResultVector.emplace_back(bufferString);
	} // if
} // function
#endif

/* Auxiliary function for printing a m128 value
 */
void vDump__m128i( __m128i value )
{
	union
	{
		__m128i value;
		int16_t arr[8];
	} converter;
	converter.value = value;

	for ( int iterator = 0; iterator < 8; iterator++ )
	{
		std::cout << converter.arr[iterator] << "\t";
	} // for

	std::cout << std::endl;
} // function

/* Log2 computation for int32 values.
 */
static inline int int_log2( uint32_t v )
{
	int c = 0;
	if ( v & 0xffff0000u ) { v >>= 16; c |= 16; }
	if ( v & 0xff00 ) { v >>= 8; c |= 8; }
	if ( v & 0xf0 ) { v >>= 4; c |= 4; }
	if ( v & 0xc ) { v >>= 2; c |= 2; }
	if ( v & 0x2 ) c |= 1;
	return c;
} // function

/* Little helper that dumps a string to a file using C++ strings.
 */
void vWriteStringToFile( const char* pcFileName, std::string &sString )
{
	std::ofstream xFileOutputStream;
	xFileOutputStream.open( pcFileName );
	xFileOutputStream << sString;
	xFileOutputStream.close();
} // function

/* Little helper function for showing the first lines of some huge file.
 */
void cShowFirstLinesOfFile()
{
	std::string sString;
	std::ifstream infile;
	infile.open( "" );
	for ( size_t uiCounter = 0; uiCounter < 10 && !infile.eof(); uiCounter++ ) // To get you all the lines.
	{
		getline( infile, sString ); // Saves the line in STRING.
		std::cout << sString << "\n"; // Prints our STRING.
	} // while
	infile.close();
} // function

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

/* Heat diagram color computation scheme.
 */
RGB_Color GetColour( double dValue, double dMinValue, double dMaxValue )
{
	RGB_Color xRGB_Color( 1.0, 1.0, 1.0 ); // white
	double dRange;

	if ( dValue < dMinValue )
		dValue = dMinValue;
	if ( dValue > dMaxValue )
		dValue = dMaxValue;
	
	dRange = dMaxValue - dMinValue;

	if ( dValue < (dMinValue + 0.25 * dRange) )
	{
		xRGB_Color.dRed = 0;
		xRGB_Color.dGreen = 4 * (dValue - dMinValue) / dRange;
	} // if
	else if ( dValue < (dMinValue + 0.5 * dRange) )
	{
		xRGB_Color.dRed = 0;
		xRGB_Color.dBlue = 1 + 4 * (dMinValue + 0.25 * dRange - dValue) / dRange;
	} // else if
	else if ( dValue < (dMinValue + 0.75 * dRange) )
	{
		xRGB_Color.dRed = 4 * (dValue - dMinValue - 0.5 * dRange) / dRange;
		xRGB_Color.dBlue = 0;
	} // else if
	else
	{
		xRGB_Color.dGreen = 1 + 4 * (dMinValue + 0.75 * dRange - dValue) / dRange;
		xRGB_Color.dBlue = 0;
	} // else

	return( xRGB_Color );
} // function


/* ******************************************
 * Deprecated C I/O functions originating from BWA code
 * ******************************************
 */
#include <stdarg.h> // for GCC

void err_fatal( const char *header, const char *fmt, ... )
{
	va_list args;
	va_start( args, fmt );
	fprintf( stderr, "[%s] ", header );
	vfprintf( stderr, fmt, args );
	fprintf( stderr, "\n" );
	va_end( args );
	exit( EXIT_FAILURE );
}

void err_fatal_core( const char *header, const char *fmt, ... )
{
	va_list args;
	va_start( args, fmt );
	fprintf( stderr, "[%s] ", header );
	vfprintf( stderr, fmt, args );
	fprintf( stderr, " Abort!\n" );
	va_end( args );
	abort();
}

void _err_fatal_simple( const char *func, const char *msg )
{
	fprintf( stderr, "[%s] %s\n", func, msg );
	exit( EXIT_FAILURE );
}

void _err_fatal_simple_core( const char *func, const char *msg )
{
	fprintf( stderr, "[%s] %s Abort!\n", func, msg );
	abort();
}

FILE *err_xopen_core(const char *func, const char *fn, const char *mode)
{
	FILE *fp = 0;
	if (strcmp(fn, "-") == 0)
		return (strstr(mode, "r"))? stdin : stdout;
	if ((fp = fopen(fn, mode)) == 0) {
		err_fatal(func, "fail to open file '%s' : %s", fn, strerror(errno));
	}
	return fp;
}

FILE *err_xreopen_core(const char *func, const char *fn, const char *mode, FILE *fp)
{
	if (freopen(fn, mode, fp) == 0) {
		err_fatal(func, "fail to open file '%s' : %s", fn, strerror(errno));
	}
	return fp;
}

gzFile err_xzopen_core(const char *func, const char *fn, const char *mode)
{
	gzFile fp;
	if (strcmp(fn, "-") == 0) {
#ifdef _MSC_VER
		fp = gzdopen(_fileno((strstr(mode, "r"))? stdin : stdout), mode);
#else
		fp = gzdopen(fileno((strstr(mode, "r"))? stdin : stdout), mode);
#endif
		/* According to zlib.h, this is the only reason gzdopen can fail */
		if (!fp) err_fatal(func, "Out of memory");
		return fp;
	}
	/* gzopen : if the file is not compressed, gzopen opens without compression.
	 */
	if ((fp = gzopen(fn, mode)) == 0) {
		err_fatal(func, "fail to open file '%s' : %s", fn, errno ? strerror(errno) : "Out of memory");
	}
	return fp;
}

size_t err_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fwrite(ptr, size, nmemb, stream);
	if (ret != nmemb) 
		_err_fatal_simple("fwrite", strerror(errno));
	return ret;
}

size_t err_fread_noeof(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fread(ptr, size, nmemb, stream);
	if (ret != nmemb)
	{
		_err_fatal_simple("fread", ferror(stream) ? strerror(errno) : "Unexpected end of file");
	}
	return ret;
}

int err_gzread(gzFile file, void *ptr, unsigned int len)
{
	int ret = gzread(file, ptr, len);

	if (ret < 0)
	{
		int errnum = 0;
		const char *msg = gzerror(file, &errnum);
		_err_fatal_simple("gzread", Z_ERRNO == errnum ? strerror(errno) : msg);
	}

	return ret;
}

int err_fseek(FILE *stream, long offset, int whence)
{
	int ret = fseek(stream, offset, whence);
	if (0 != ret)
	{
		_err_fatal_simple("fseek", strerror(errno));
	}
	return ret;
}

long err_ftell(FILE *stream)
{
	long ret = ftell(stream);
	if (-1 == ret)
	{
		_err_fatal_simple("ftell", strerror(errno));
	}
	return ret;
}

int err_printf(const char *format, ...) 
{
	va_list arg;
	int done;
	va_start(arg, format);
	done = vfprintf(stdout, format, arg);
	int saveErrno = errno;
	va_end(arg);
	if (done < 0) _err_fatal_simple("vfprintf(stdout)", strerror(saveErrno));
	return done;
}

int err_fprintf(FILE *stream, const char *format, ...) 
{
	va_list arg;
	int done;
	va_start(arg, format);
	done = vfprintf(stream, format, arg);
	int saveErrno = errno;
	va_end(arg);
	if (done < 0) _err_fatal_simple("vfprintf", strerror(saveErrno));
	return done;
}

int err_fputc(int c, FILE *stream)
{
	int ret = putc(c, stream);
	if (EOF == ret)
	{
		_err_fatal_simple("fputc", strerror(errno));
	}

	return ret;
}

int err_fputs(const char *s, FILE *stream)
{
	int ret = fputs(s, stream);
	if (EOF == ret)
	{
		_err_fatal_simple("fputs", strerror(errno));
	}

	return ret;
}

int err_fflush(FILE *stream) 
{
    int ret = fflush(stream);
    if (ret != 0) _err_fatal_simple("fflush", strerror(errno));

#ifdef FSYNC_ON_FLUSH
	/* Calling fflush() ensures that all the data has made it to the
	   kernel buffers, but this may not be sufficient for remote filesystems
	   (e.g. NFS, lustre) as an error may still occur while the kernel
	   is copying the buffered data to the file server.  To be sure of
	   catching these errors, we need to call fsync() on the file
	   descriptor, but only if it is a regular file.  */
	{
		struct stat sbuf;
		if (0 != fstat(fileno(stream), &sbuf))
			_err_fatal_simple("fstat", strerror(errno));
		
		if (S_ISREG(sbuf.st_mode))
		{
			if (0 != fsync(fileno(stream)))
				_err_fatal_simple("fsync", strerror(errno));
		}
	}
#endif
    return ret;
}

int err_fclose(FILE *stream) 
{
	int ret = fclose(stream);
	if (ret != 0) _err_fatal_simple("fclose", strerror(errno));
	return ret;
}

int err_gzclose(gzFile file)
{
	int ret = gzclose(file);
	if (Z_OK != ret)
	{
		_err_fatal_simple("gzclose", Z_ERRNO == ret ? strerror(errno) : zError(ret));
	}

	return ret;
}