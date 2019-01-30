/**
 * @file support.cpp
 * @author Arne Kutzner
 */
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#include <direct.h>
#endif


#include <cerrno>
#include <cstring>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

#include "util/exception.h"
#include "util/support.h"
/* Constructs the full file name for some prefix, suffix combination.
 * Returns by value for convenience purposes.
 */
std::string fullFileName( const char* pcFileNamePrefix, const char* pcSuffix )
{
    std::string sFileName( pcFileNamePrefix );
    sFileName.push_back( '.' );
    return sFileName.append( pcSuffix );
} // method

bool EXPORTED fileExists( const std::string& rsFile )
{
    struct stat buffer;
    return ( stat( rsFile.c_str( ), &buffer ) == 0 );
} // function

bool EXPORTED ends_with( const std::string& rsX, const std::string& rsEnd )
{
    return rsX.compare( rsX.length( ) - rsEnd.length( ), rsEnd.length( ), rsEnd ) == 0;
} // function

void makeDir( const std::string& rsFile )
{
    int nError = 0;
#if defined( _WIN32 )
    nError = _mkdir( rsFile.c_str( ) ); // can be used on Windows
#else
    mode_t nMode = 0733; // UNIX style permissions
    nError = mkdir( rsFile.c_str( ), nMode ); // can be used on non-Windows
#endif
    if( nError != 0 )
    {
        // 17 == file exists we want to ignore this error
        if( errno != 17 )
        {
            throw AnnotatedException( std::string( "Could not create Dir: " )
                                          .append( rsFile )
                                          .append( " errno: " )
                                          .append( std::strerror( errno ) )
                                          .c_str( ) );
        } // if
    } // if
} // function

// taken from: https://stackoverflow.com/questions/281818/unmangling-the-result-of-stdtype-infoname
#ifdef __GNUG__
#include <cstdlib>
#include <cxxabi.h>
#include <memory>
std::string EXPORTED demangle( const char* name )
{
    int status = -4; // some arbitrary value to eliminate the compiler warning
    std::unique_ptr<char, void ( * )( void* )> res{abi::__cxa_demangle( name, NULL, NULL, &status ), std::free};
    return ( status == 0 ) ? res.get( ) : name;
} // function
#else

std::string EXPORTED demangle( const char* name )
{
    return std::string( name );
} // function

#endif