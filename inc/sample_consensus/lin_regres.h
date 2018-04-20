#include <iostream>
#include <vector>
#include <cassert>

template<typename TP>
TP mean( const std::vector<TP> &vVector )
{
	TP sum = 0;
	for (size_t i = 0; i < vVector.size(); i++)
		sum = sum + vVector[i];
	return sum / (TP)vVector.size();
} // function

/* Corrected sum of squares
 */
template<typename TP>
TP corrSumOfSquares( const std::vector<TP> &vaVector,
				     std::vector<TP> &vdVector,
					 const TP mean )
{	
	assert( vaVector.size() == vdVector.size() );

	TP sum = 0;
	for (size_t i = 0; i < vaVector.size(); i++) {
		vdVector[i] = vaVector[i] - mean;
		sum = sum + (vdVector[i] * vdVector[i]);
	}
	return sum;
} // function


template<typename TP>
std::pair<TP, TP> lin_regres( const std::vector<TP> &vxVec, 
							  const std::vector<TP> &vyVec )
{
	assert( vxVec.size() == vyVec.size() );
	
	size_t uiNumSamples = vxVec.size();

	std::vector<TP> vdxVector( uiNumSamples ); // expensive
	std::vector<TP> vdyVector( uiNumSamples ); // expensive
	
	/* mean value x: x-bar */
	auto mean_x = mean( vxVec );
	//// std::cout << "mean X: " << mean_x << std::endl;
	/* mean value x: y-bar */
	auto mean_y = mean( vyVec );
	//// std::cout << "mean Y: " << mean_y << std::endl;

	/* Corrected sum of squares S_{xx} */
	auto sx = corrSumOfSquares( vxVec, vdxVector, mean_x );
	//// std::cout << "corrSumOfSquares X: " << sx << std::endl;
	/* Corrected sum of squares S_{yy} */
	auto sy = corrSumOfSquares( vyVec, vdyVector, mean_y );
	//// std::cout << "corrSumOfSquaresn Y: " << sy << std::endl;

	/* Corrected sum of cross products S_{xy}.
	 * Can be computed without auxiliary vectors.
	 */
	TP sum_xy = 0;
	for (size_t i = 0; i < uiNumSamples; i++)
		sum_xy = sum_xy + vdxVector[i] * vdyVector[i];
	
	/* Not used for our purpose */
	//// TP deviation_x = sqrt( sx ) / uiNumSamples; 
	//// TP deviation_y = sqrt( sy ) / uiNumSamples;
	//// TP corr_coff = sum_xy / (uiNumSamples * deviation_x * deviation_y);
	//// TP reg_coff_xy = corr_coff * (deviation_x / deviation_y);
	//// TP reg_coff_yx = corr_coff * (deviation_y / deviation_x);
	
	/* slope (b): b = S_{ xy } / S_{xx} */
	TP slope = sum_xy / sx;
	//// std::cout << "slope: " << slope << std::endl;

	/* intercept (a): a = mean(Y) - b * mean(X) */
	TP intercept = mean_y - slope * mean_x;
	//// std::cout << "intercept: " << intercept << std::endl;
	
	return std::pair<TP, TP>(slope, intercept);
} // function

std::pair<double, double> lin_regres_double( const std::vector<double> &vA,
											 const std::vector<double> &vB )
{
	return lin_regres( vA, vB );
} // function

void test_lin_regres()
{
	std::vector<double> vX = { 1, 2, 3, 4, 5 };
	std::vector<double> vY = { 5, 7, 9, 6, 8 };

	auto xSlopeIntercept = lin_regres_double( vX, vY );

	std::cout << "slope: " << xSlopeIntercept.first << std::endl;
	std::cout << "intercept: " << xSlopeIntercept.second << std::endl;
}