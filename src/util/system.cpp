#include "util/system.h"


// taken from: https://stackoverflow.com/questions/281818/unmangling-the-result-of-stdtype-infoname
#ifdef __GNUG__
#include <cstdlib>
#include <cxxabi.h>
#include <memory>
std::string demangle( const char* name )
{
    int status = -4; // some arbitrary value to eliminate the compiler warning
    std::unique_ptr<char, void ( * )( void* )> res{abi::__cxa_demangle( name, NULL, NULL, &status ), std::free};
    return ( status == 0 ) ? res.get( ) : name;
} // function
#else
// does nothing if not g++
std::string demangle( const char* name )
{
    return name;
} // function
#endif