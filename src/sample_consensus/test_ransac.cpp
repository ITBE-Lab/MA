#include "sample_consensus/test_ransac.h"
#include "module/module.h"
#include "sample_consensus/lin_regres.h"
#include "sample_consensus/ransac.h"
#include "sample_consensus/sac_model_line.h"

using namespace libMA;

std::pair<double, double> run_ransac( const std::vector<double>& rvxValues,
                                      const std::vector<double>& rvyValues,
                                      /*std::shared_ptr<Seeds> pvIngroup,*/
                                      double fDMA )
{
    /* Fill the PointCloud with the vector data
     */
    sensor_msgs::PointCloud xCloud;
    for ( size_t uiItr = 0; uiItr < rvxValues.size( ); uiItr++ )
    {
        xCloud.points.push_back(
            geometry_msgs::Point32{rvxValues[ uiItr ], rvyValues[ uiItr ], 0} );
    } // for
    // std::cout << "Process " << xCloud.points.size() << " points" << std::endl;

    /* Create the xModel and set the DataSet fot the xModel
     */
    sample_consensus::SACModelLine xModel;
    xModel.setDataSet( &xCloud );

    // std::cout << "cpp residual_threshold is: " << dMAD << std::endl;

    /* dMAD argument is the threshold used for identifying inliers.
     * t – threshold value to determine when a data point fits a model
     * In the SciPy Documentation:
     *    residual_threshold : float, optional
     *    Maximum residual for a data sample to be classified as an inlier.
     *    By default the threshold is chosen as the MAD (median absolute deviation) of the target
     * values y.
     */
    sample_consensus::RANSAC xRansac( &xModel, fDMA );

    /* Do the inliers extraction by performing an iterative process.
     * This line is equal to: ransac.fit(X, y) in the Pyhton code.
     * Argument is the debug level.
     */
    bool bSuccessfulModel = xRansac.computeModel( 0 ); // debug level 0

    if ( !bSuccessfulModel )
    {
        DEBUG( std::cerr << "ransac failed to find ingroup" << std::endl; ) // DEBUG
        // return that no slope could be determined
        return std::make_pair( std::nan( "" ), std::nan( "" ) );
    } // if

    //  /* Print the coefficient ids (of 2 anchor points) */
    //  for(auto fX : xModel.getBestModel())
    //  	std::cout << fX << std::endl;
    //
    //  /* Print the coefficient data (of 2 anchor points) */
    //  for(auto fX : xModel.getModelCoefficients())
    //      std::cout << fX << std::endl;
    //
    //  /* Print the set of all dots that remain as good points */
    //  for(auto fX : xModel.getBestInliers())
    //  {
    //  	std::cout << "(" << rvxValues[fX] << "," << rvyValues[fX] << ")";
    //  }// for
    //  std::cout << std::endl;

    /*
     * We have now a model with 900 inliers.
     * These inliers have to be used for the computation of the linear regression.
     * There seems to be no code in the ransac stuff that performs linear regression.
     */

    // std::cout << "NOW WE HAVE TO DO LINEAR REGRESSION..." << std::endl;

    std::vector<double> vX;
    vX.reserve( xModel.getBestInliers( ).size( ) );
    std::vector<double> vY;
    vY.reserve( xModel.getBestInliers( ).size( ) );

    for ( int iIndex : xModel.getBestInliers( ) )
    {
        vX.push_back( xCloud.points[ iIndex ].x );
        vY.push_back( xCloud.points[ iIndex ].y );
        // DEBUG(
        //    pvIngroup->push_back(Seed(xCloud.points[iIndex].y, 1, xCloud.points[iIndex].x));
        //)
    } // for

    auto xSlopeIntercept = lin_regres_double( vX, vY );

    // std::cout << "slope: " << xSlopeIntercept.first << std::endl;
    // std::cout << "intercept: " << xSlopeIntercept.second << std::endl;

    return xSlopeIntercept;
} // function

void test_ransac( )
{
    std::vector<double> vX = {1, 2, 3, 4, 5};
    std::vector<double> vY = {5, 7, 9, 6, 8};
    double fMAD = medianAbsoluteDeviation<double>( vY );

    auto xSlopeIntercept = run_ransac( vX, vY, /*std::make_shared<Seeds>(),*/ fMAD );

    std::cout << "slope: " << xSlopeIntercept.first << std::endl;
    std::cout << "intercept: " << xSlopeIntercept.second << std::endl;
} // function

#ifdef WITH_PYTHON
void export_ransac( )
{
    boost::python::def( "test_ransac", &test_ransac );
} // function
#endif
