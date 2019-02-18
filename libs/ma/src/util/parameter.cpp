#include "util/parameter.h"


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

// provide implementation for forward declared constructor.
EXPORTED AlignerParameterBase::AlignerParameterBase( Presetting* pPresetting, const std::string& sName,
                                            const std::string& sDescription )
    : sName( sName ), sDescription( sDescription )
{
    if( pPresetting != nullptr )
        pPresetting->registerParameter( this );
} // constructor