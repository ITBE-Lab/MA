#include <cassert>
#include <iostream>
#include <vector>

#define DO_ERROR ( 0 )

template <typename TP> TP mean( const std::vector<TP> &vVector )
{
    TP sum = 0;
    for( size_t i = 0; i < vVector.size( ); i++ )
        sum = sum + vVector[ i ];
    return sum / (TP)vVector.size( );
} // function

/* Corrected sum of squares
 * FUNCTION HAS SIDE EFFECT
 */
template <typename TP>
TP corrSumOfSquares( const std::vector<TP> &vaVector, std::vector<TP> &vdVector, const TP mean )
{
    assert( vaVector.size( ) == vdVector.size( ) );

    TP sum = 0;
    for( size_t i = 0; i < vaVector.size( ); i++ )
    {
        vdVector[ i ] = vaVector[ i ] - mean;
        sum = sum + ( vdVector[ i ] * vdVector[ i ] );
    }
    return sum;
} // function

// CAREFULL HAS SIDEFFECTS
template <typename TP>
TP variance( const std::vector<TP> &vaVector, std::vector<TP> &vdVector, const TP mean )
{
    return corrSumOfSquares( vaVector, vdVector, mean ) / vaVector.size( );
}

// CAREFULL HAS SIDEFFECTS
template <typename TP>
TP covariance( const std::vector<TP> &vaVector, const std::vector<TP> &vdVector, const TP mean_x,
               const TP mean_y )
{
    double sum = 0;
    for( size_t i = 0; i < vaVector.size( ); i++ )
    {
        sum += ( vaVector[ i ] - mean_x ) * ( vdVector[ i ] - mean_y );
    } // for
    return sum / vaVector.size( );
} // function

template <typename TP>
std::pair<TP, TP> lin_regres( const std::vector<TP> &vxVec, const std::vector<TP> &vyVec )
{
    assert( vxVec.size( ) == vyVec.size( ) );

    size_t uiNumSamples = vxVec.size( );

    std::vector<TP> vdxVector( uiNumSamples ); // expensive
    std::vector<TP> vdyVector( uiNumSamples ); // expensive

    /* mean value x: x-bar */
    auto mean_x = mean( vxVec );
    // std::cout << "mean X: " << mean_x << std::endl;
    /* mean value x: y-bar */
    auto mean_y = mean( vyVec );
    // std::cout << "mean Y: " << mean_y << std::endl;

    /* Corrected sum of squares S_{xx} */
    auto sx = corrSumOfSquares( vxVec, vdxVector, mean_x );
    //// std::cout << "corrSumOfSquares X: " << sx << std::endl;
    /* Corrected sum of squares S_{yy} */
#if DO_ERROR
    auto sy =
#endif
        /*auto sy =*/corrSumOfSquares( vyVec, vdyVector, mean_y );
    //// std::cout << "corrSumOfSquaresn Y: " << sy << std::endl;

#if DO_ERROR
    auto varX = sx / vdxVector.size( );
    auto varY = sy / vdyVector.size( );

    auto cov = covariance( vxVec, vyVec, mean_x, mean_y );
#endif

    /* Corrected sum of cross products S_{xy}.
     * Can be computed without auxiliary vectors.
     */
    TP sum_xy = 0;
    for( size_t i = 0; i < uiNumSamples; i++ )
        sum_xy = sum_xy + vdxVector[ i ] * vdyVector[ i ];

    /* Not used for our purpose */
    //// TP deviation_x = sqrt( sx ) / uiNumSamples;
    //// TP deviation_y = sqrt( sy ) / uiNumSamples;
    //// TP corr_coff = sum_xy / (uiNumSamples * deviation_x * deviation_y);
    //// TP reg_coff_xy = corr_coff * (deviation_x / deviation_y);
    //// TP reg_coff_yx = corr_coff * (deviation_y / deviation_x);

    /* slope (b): b = S_{ xy } / S_{xx} */
    TP slope = sum_xy / sx;
    // std::cout << "slope: " << slope << std::endl;

    /* intercept (a): a = mean(Y) - b * mean(X) */
    TP intercept = mean_y - slope * mean_x;
    // std::cout << "intercept: " << intercept << std::endl;

    /// TP reg_ss = (sum_xy*sum_xy) / sx;

    TP total_ss = 0;
    for( auto y : vyVec )
        total_ss += ( y - mean_y ) * ( y - mean_y );

        // old error
        //// TP error = ( total_ss - reg_ss ) / (vyVec.size() - 2);
#if DO_ERROR
    double r_den = std::sqrt( varX * varY );
    double r = 0;
    if( r_den == 0.0 )
        r = 0.0;
    else
    {
        r = cov / r_den;
        // test for numerical error propagation
        if( r > 1.0 )
            r = 1.0;
        else if( r < -1.0 )
            r = -1.0;
    } // else
    TP error = std::sqrt( ( 1 - r * r ) * varY / varX / ( vxVec.size( ) - 2 ) );
    if( vyVec.size( ) > 2 )
        std::cout << "lin reg ERROR: " << error << std::endl;
#endif

    return std::pair<TP, TP>( std::atan( slope ), -intercept / slope );
} // function

std::pair<double, double> lin_regres_double( const std::vector<double> &vA,
                                             const std::vector<double> &vB )
{
    return lin_regres( vA, vB );
} // function

void test_lin_regres( )
{
    std::vector<double> vX = {1, 2, 3, 4, 5};
    std::vector<double> vY = {5, 7, 9, 6, 8};

    auto xSlopeIntercept = lin_regres_double( vX, vY );

    std::cout << "slope: " << xSlopeIntercept.first << std::endl;
    std::cout << "intercept: " << xSlopeIntercept.second << std::endl;
}