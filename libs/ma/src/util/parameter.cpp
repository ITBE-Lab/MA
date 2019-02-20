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

template <> EXPORTED bool genericStringToValue<bool>( const std::string& sString )
{
    return sString == "true";
} // function


const EXPORTED std::pair<size_t, std::string> Presetting::DP_PARAMETERS = std::make_pair( 5, "Dynamic Programming");
const EXPORTED std::pair<size_t, std::string> Presetting::HEURISTIC_PARAMETERS = std::make_pair( 4, "Heuristics");
const EXPORTED std::pair<size_t, std::string> Presetting::SEEDING_PARAMETERS = std::make_pair( 1, "Seeding");
const EXPORTED std::pair<size_t, std::string> Presetting::SOC_PARAMETERS = std::make_pair( 2, "Strip of Consideration");
const EXPORTED std::pair<size_t, std::string> Presetting::PAIRED_PARAMETERS = std::make_pair( 0, "Paired Reads");
const EXPORTED std::pair<size_t, std::string> Presetting::SAM_PARAMETERS = std::make_pair( 3, "SAM Output");