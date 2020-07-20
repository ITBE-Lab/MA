/**
 * @file abstractFilter.h
 * @author Markus Schmidt
 */

#pragma once

namespace libMSV
{

#define ANALYZE_FILTERS ( 1 )

class AbstractFilter
{
  public:
#if ANALYZE_FILTERS
    std::string sName;
    size_t uiFilterKept = 0;
    size_t uiFilterTotal = 0;
    std::mutex xLock;
    static bool bSilent;
#endif
    AbstractFilter( std::string sName )
#if ANALYZE_FILTERS
        : sName( sName )
#endif
    {} // constructor
#if ANALYZE_FILTERS

    ~AbstractFilter( )
    {
        if( !bSilent && uiFilterTotal > 0 )
            std::cout << std::fixed << std::setprecision( 2 ) << "~" << sName << ": filter kept and eliminated "
                      << uiFilterKept << " and " << uiFilterTotal - uiFilterKept
                      << " elements respectiveley.\n\tThat's " << 100.0 * uiFilterKept / (double)uiFilterTotal
                      << "% and " << 100.0 - 100.0 * uiFilterKept / (double)uiFilterTotal << "% respectiveley."
                      << std::endl
                      << std::defaultfloat;
    } // deconstructor
#endif
}; // class

} // namespace libMSV