#include "util/parameter.h"
#include <algorithm>

template <> EXPORTED std::string genericStringToValue<std::string>( const std::string& sString )
{
    return std::string( sString );
} // function

template <> EXPORTED int genericStringToValue<int>( const std::string& sString )
{
    return stoi( sString );
} // function

template <> EXPORTED bool genericStringToValue<bool>( const std::string& sString )
{
    return sString == "true";
} // function