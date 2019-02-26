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
    return sString == "true";
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


const EXPORTED std::pair<size_t, std::string> Presetting::DP_PARAMETERS = std::make_pair( 5, "Dynamic Programming" );
const EXPORTED std::pair<size_t, std::string> Presetting::HEURISTIC_PARAMETERS = std::make_pair( 4, "Heuristics" );
const EXPORTED std::pair<size_t, std::string> Presetting::SEEDING_PARAMETERS = std::make_pair( 1, "Seeding" );
const EXPORTED std::pair<size_t, std::string> Presetting::SOC_PARAMETERS =
    std::make_pair( 2, "Strip of Consideration" );
const EXPORTED std::pair<size_t, std::string> Presetting::PAIRED_PARAMETERS = std::make_pair( 0, "Paired Reads" );
const EXPORTED std::pair<size_t, std::string> Presetting::SAM_PARAMETERS = std::make_pair( 3, "SAM Output" );
const EXPORTED std::pair<size_t, std::string> Presetting::SV_PARAMETERS =
    std::make_pair( 6, "Structural Variants Caller" );
const EXPORTED std::pair<size_t, std::string> GeneralParameter::GENERAL_PARAMETER =
    std::make_pair( 0, "General Parameter" );