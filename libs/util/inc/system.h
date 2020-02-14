/**
 * @file system.h
 * @brief Implements some analysis functions.
 * @author Markus Schmidt
 */
#pragma once

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <chrono> // required for getting the high resolution clock
#include <iostream>
#include <string>

#if __GNUC__
#include <time.h>
#endif
/// @endcond

#define defDO_LOG true
#define defSUPRESS_LOG false

/* Meta function for measuring durations of function executions.
 * Using count() applied to the returned object we get the time in seconds.
 */
template <class FUNCTOR> std::chrono::duration<double> metaMeasureDuration( FUNCTOR&& f )
{ /* record start time
   */
    auto start = std::chrono::high_resolution_clock::now( );
    f( );

    /* record end time
     */
    auto end = std::chrono::high_resolution_clock::now( );
    return end - start;
} // meta function

template <bool bLog, class FUNCTOR>
double metaMeasureAndLogDuration( const std::string& sLogText, // additional logging text
                                  FUNCTOR&& f // the functor called for measuring execution time
)
{
    if( bLog )
    {
        // Measure duration and log.
        auto xDuration = metaMeasureDuration( std::forward<FUNCTOR>( f ) );
        std::cout << sLogText << " required " << xDuration.count( ) * 1000 << " milliseconds." << std::endl;
        return xDuration.count( ) * 1000;
    } // if
    else
    {
        // Simply call the functor.
        f( );
        return 0;
    } // else
} // function

#if _MSC_VER
template <class TP_FUNC_APPLY> __int64 time_call_( TP_FUNC_APPLY&& f )
{
    __int64 begin = GetTickCount( );
    f( );
    return GetTickCount( ) - begin;
} // function

#elif __GNUC__
template <class TP_FUNC_APPLY> timespec time_call_( TP_FUNC_APPLY&& f )
{
    timespec startTime, endTime, differenceTime;

    /* the function call embedded by two time measurements
     */
    clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &startTime );
    f( );
    clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &endTime );

    /* compute the time difference
     */
    if( ( endTime.tv_nsec - startTime.tv_nsec ) < 0 )
    {
        /* overflow occurred
         */
        differenceTime.tv_sec = endTime.tv_sec - startTime.tv_sec - 1;
        differenceTime.tv_nsec = 1000000000 + endTime.tv_nsec - startTime.tv_nsec;
    }
    else
    {
        differenceTime.tv_sec = endTime.tv_sec - startTime.tv_sec;
        differenceTime.tv_nsec = endTime.tv_nsec - startTime.tv_nsec;
    }

    return differenceTime;
} // function
#endif

#ifdef __GNUC__
std::string timeAsString( const timespec& time );
#endif

#ifdef _MSC_VER
std::string timeAsString( const __int64& time );
#endif
