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

template <> EXPORTED double genericStringToValue<double>( const std::string& sString )
{
    return stod( sString );
} // function

template <> EXPORTED float genericStringToValue<float>( const std::string& sString )
{
    return (float)stod( sString );
} // function

template <> EXPORTED bool genericStringToValue<bool>( const std::string& sString )
{
    std::string sTmp = sString;
    std::locale loc;
    for( std::string::size_type i = 0; i < sTmp.length( ); ++i )
        sTmp[ i ] = std::tolower( sTmp[ i ], loc );
    if( sTmp == "true" )
        return true;
    if( sTmp == "false" )
        return true;
    throw std::runtime_error( "Boolean flags accept no other values than 'true' and 'false'." );
} // function

template <> EXPORTED std::string AlignerParameter<bool>::asText( ) const
{
    if( this->get( ) )
        return "true";
    return "false";
} // method

template <> EXPORTED std::string AlignerParameter<double>::asText( ) const
{
    // Taken from: https://stackoverflow.com/questions/13686482/c11-stdto-stringdouble-no-trailing-zeros
    std::string sText( std::to_string( value ) );
    int iOffset = 1;
    if( sText.find_last_not_of( '0' ) == sText.find( '.' ) )
        iOffset = 0;
    sText.erase( sText.find_last_not_of( '0' ) + iOffset, std::string::npos );
    return sText;
} // method

template <> EXPORTED std::string AlignerParameter<float>::asText( ) const
{
    // Taken from: https://stackoverflow.com/questions/13686482/c11-stdto-stringdouble-no-trailing-zeros
    std::string sText( std::to_string( value ) );
    int iOffset = 1;
    if( sText.find_last_not_of( '0' ) == sText.find( '.' ) )
        iOffset = 0;
    sText.erase( sText.find_last_not_of( '0' ) + iOffset, std::string::npos );
    return sText;
} // method

template <> std::string AlignerParameter<int>::asText( ) const
{
    return std::to_string( this->get( ) );
} // method

template <> std::string AlignerParameter<std::string>::asText( ) const
{
    return this->get( );
} // method
