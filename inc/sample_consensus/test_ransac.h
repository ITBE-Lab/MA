
#ifndef TEST_RANSAC_H
#define TEST_RANSAC_H

#include "container/seed.h"
#include <algorithm>

/* Function to find mean of the array elements.
 */
template <typename TP> TP Mean( const std::vector<TP> arr )
{
    // Calculate sum of all elements.
    TP sum = 0;
    for( size_t i = 0; i < arr.size( ); i++ )
        sum += arr[ i ];
    return sum / arr.size( );
} // function

/* Function to find mean of the array elements.
 */
template <typename TP> TP Median( std::vector<TP> arr )
{
    std::sort( arr.begin( ), arr.end( ) );
    // Calculate sum of all elements.
    if( arr.size( ) == 0 )
    {
        return 0; // Undefined, really.
    }
    else if( arr.size( ) == 1 )
    {
        return arr[ 0 ];
    }

    if( arr.size( ) % 2 == 0 )
    {
        return ( arr[ arr.size( ) / 2 - 1 ] + arr[ arr.size( ) / 2 ] ) / 2;
    }
    return arr[ arr.size( ) / 2 ];
} // function

// Function to find mean absolute
// deviation of given elements.
template <typename TP> TP meanAbsoluteDeviation( const std::vector<TP> arr )
{
    // Calculate the sum of absolute
    // deviation about mean.
    TP absSum = 0;
    TP mean = Mean( arr );
    for( size_t i = 0; i < arr.size( ); i++ )
        absSum += abs( arr[ i ] - mean );

    // Return mean absolute deviation about mean.
    return absSum / arr.size( );
} // function

// Function to find mean absolute
// deviation of given elements.
template <typename TP> TP medianAbsoluteDeviation( std::vector<TP> arr )
{
    // Calculate the sum of absolute
    // deviation about mean.
    // Python: np.median(np.abs(y - np.median(y)))


    TP median = Median( arr );
    std::vector<TP> vDevArr;
    for( size_t i = 0; i < arr.size( ); i++ )
        if( arr[ i ] - median < 0 )
            vDevArr.push_back( -( arr[ i ] - median ) );
        else
            vDevArr.push_back( arr[ i ] - median );

    // Return mean absolute deviation about mean.
    return Median( vDevArr );
} // function

std::pair<double, double> run_ransac( const std::vector<double> &rvxValues,
                                      const std::vector<double> &rvyValues,
                                      /*std::shared_ptr<libMA::Seeds> pvIngroup,*/
                                      double fDMA );

void export_ransac( );

#endif