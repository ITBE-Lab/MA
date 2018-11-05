/*
 * Copyright (c) 2008 Radu Bogdan Rusu <rusu -=- cs.tum.edu>
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * $Id: sac_model_plane.cpp 21050 2009-08-07 21:24:30Z jfaustwg $
 *
 */

/** \author Radu Bogdan Rusu */

#include <sample_consensus/sac_model_plane.h>
// #include <point_cloud_mapping/geometry/nearest.h>

namespace sample_consensus
{
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Get 3 random non-collinear points as data samples and return them as point indices.
 * \param iterations the internal number of iterations used by SAC methods
 * \param samples the resultant model samples
 * \note assumes unique points!
 */
void SACModelPlane::getSamples( int &iterations, std::vector<int> &samples )
{
    samples.resize( 3 );
    double trand = indices_.size( ) / ( RAND_MAX + 1.0 );

    // Get a random number between 1 and max_indices
    int idx = (int)( rand( ) * trand );
    // Get the index
    samples[ 0 ] = indices_.at( idx );

    // Get a second point which is different than the first
    do
    {
        idx = (int)( rand( ) * trand );
        samples[ 1 ] = indices_.at( idx );
        iterations++;
    } while( samples[ 1 ] == samples[ 0 ] );
    iterations--;

    double Dx1, Dy1, Dz1, Dx2, Dy2, Dz2, Dy1Dy2;
    // Compute the segment values (in 3d) between XY
    Dx1 = cloud_->points[ samples[ 1 ] ].x - cloud_->points[ samples[ 0 ] ].x;
    Dy1 = cloud_->points[ samples[ 1 ] ].y - cloud_->points[ samples[ 0 ] ].y;
    Dz1 = cloud_->points[ samples[ 1 ] ].z - cloud_->points[ samples[ 0 ] ].z;

    int iter = 0;
    do
    {
        // Get the third point, different from the first two
        do
        {
            idx = (int)( rand( ) * trand );
            samples[ 2 ] = indices_.at( idx );
            iterations++;
        } while( ( samples[ 2 ] == samples[ 1 ] ) || ( samples[ 2 ] == samples[ 0 ] ) );
        iterations--;

        // Compute the segment values (in 3d) between XZ
        Dx2 = cloud_->points[ samples[ 2 ] ].x - cloud_->points[ samples[ 0 ] ].x;
        Dy2 = cloud_->points[ samples[ 2 ] ].y - cloud_->points[ samples[ 0 ] ].y;
        Dz2 = cloud_->points[ samples[ 2 ] ].z - cloud_->points[ samples[ 0 ] ].z;

        Dy1Dy2 = Dy1 / Dy2;
        iter++;

        if( iter > MAX_ITERATIONS_COLLINEAR )
        {
            std::cout << "[SACModelPlane::getSamples] WARNING: Could not select 3 non collinear "
                         "points in %d iterations!"
                      << MAX_ITERATIONS_COLLINEAR << std::endl;
            break;
        }
        iterations++;
    }
    // Use Zoli's method for collinearity check
    while( ( ( Dx1 / Dx2 ) == Dy1Dy2 ) && ( Dy1Dy2 == ( Dz1 / Dz2 ) ) );
    iterations--;

    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Select all the points which respect the given model coefficients as inliers.
 * \param model_coefficients the coefficients of a plane model that we need to compute distances to
 * \param threshold a maximum admissible distance threshold for determining the inliers from the
 * outliers \param inliers the resultant model inliers \note: we should compare (e^2) < (T^2) but
 * instead we use fabs(e) because threshold is always positive \note: To get the refined inliers of
 * a model, use: ANNpoint refined_coeff = refitModel (...); selectWithinDistance (refined_coeff,
 * threshold);
 */
void SACModelPlane::selectWithinDistance( const std::vector<double> &model_coefficients,
                                          double threshold, std::vector<int> &inliers )
{
    int nr_p = 0;
    inliers.resize( indices_.size( ) );

    // Iterate through the 3d points and calculate the distances from them to the plane
    for( unsigned int i = 0; i < indices_.size( ); i++ )
    {
        // Calculate the distance from the point to the plane normal as the dot product
        // D = (P-A).N/|N|
        if( fabs( model_coefficients.at( 0 ) * cloud_->points.at( indices_.at( i ) ).x +
                  model_coefficients.at( 1 ) * cloud_->points.at( indices_.at( i ) ).y +
                  model_coefficients.at( 2 ) * cloud_->points.at( indices_.at( i ) ).z +
                  model_coefficients.at( 3 ) ) < threshold )
        {
            // Returns the indices of the points whose distances are smaller than the threshold
            inliers[ nr_p ] = indices_[ i ];
            nr_p++;
        }
    }
    inliers.resize( nr_p );
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Compute all distances from the cloud data to a given plane model.
 * \param model_coefficients the coefficients of a plane model that we need to compute distances to
 * \param distances the resultant estimated distances
 */
void SACModelPlane::getDistancesToModel( const std::vector<double> &model_coefficients,
                                         std::vector<double> &distances )
{
    distances.resize( indices_.size( ) );

    // Iterate through the 3d points and calculate the distances from them to the plane
    for( unsigned int i = 0; i < indices_.size( ); i++ )
        // Calculate the distance from the point to the plane normal as the dot product
        // D = (P-A).N/|N|
        distances[ i ] = fabs( model_coefficients.at( 0 ) * cloud_->points.at( indices_[ i ] ).x +
                               model_coefficients.at( 1 ) * cloud_->points.at( indices_[ i ] ).y +
                               model_coefficients.at( 2 ) * cloud_->points.at( indices_[ i ] ).z +
                               model_coefficients.at( 3 ) );
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Create a new point cloud with inliers projected onto the plane model.
 * \param inliers the data inliers that we want to project on the plane model
 * \param model_coefficients the *normalized* coefficients of a plane model
 * \param projected_points the resultant projected points
 */
void SACModelPlane::projectPoints( const std::vector<int> &inliers,
                                   const std::vector<double> &model_coefficients,
                                   sensor_msgs::PointCloud &projected_points )
{
    // Allocate enough space
    projected_points.points.resize( inliers.size( ) );
    projected_points.set_channels_size( (unsigned int) cloud_->get_channels_size( ) );

    // Create the channels
    for( unsigned int d = 0; d < projected_points.get_channels_size( ); d++ )
    {
        projected_points.channels[ d ].name = cloud_->channels[ d ].name;
        projected_points.channels[ d ].values.resize( inliers.size( ) );
    }

    // Get the plane normal
    // Calculate the 2-norm: norm (x) = sqrt (sum (abs (v)^2))
    // double n_norm = sqrt (model_coefficients.at (0) * model_coefficients.at (0) +
    //                      model_coefficients.at (1) * model_coefficients.at (1) +
    //                      model_coefficients.at (2) * model_coefficients.at (2));
    // model_coefficients.at (0) /= n_norm;
    // model_coefficients.at (1) /= n_norm;
    // model_coefficients.at (2) /= n_norm;
    // model_coefficients.at (3) /= n_norm;

    // Iterate through the 3d points and calculate the distances from them to the plane
    for( unsigned int i = 0; i < inliers.size( ); i++ )
    {
        // Calculate the distance from the point to the plane
        double distance_to_plane =
            model_coefficients.at( 0 ) * cloud_->points.at( inliers.at( i ) ).x +
            model_coefficients.at( 1 ) * cloud_->points.at( inliers.at( i ) ).y +
            model_coefficients.at( 2 ) * cloud_->points.at( inliers.at( i ) ).z +
            model_coefficients.at( 3 ) * 1;
        // Calculate the projection of the point on the plane
        projected_points.points[ i ].x =
            cloud_->points.at( inliers.at( i ) ).x - distance_to_plane * model_coefficients.at( 0 );
        projected_points.points[ i ].y =
            cloud_->points.at( inliers.at( i ) ).y - distance_to_plane * model_coefficients.at( 1 );
        projected_points.points[ i ].z =
            cloud_->points.at( inliers.at( i ) ).z - distance_to_plane * model_coefficients.at( 2 );
        // Copy the other attributes
        for( unsigned int d = 0; d < projected_points.get_channels_size( ); d++ )
            projected_points.channels[ d ].values[ i ] =
                cloud_->channels[ d ].values[ inliers.at( i ) ];
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Project inliers (in place) onto the given plane model.
 * \param inliers the data inliers that we want to project on the plane model
 * \param model_coefficients the *normalized* coefficients of a plane model
 */
void SACModelPlane::projectPointsInPlace( const std::vector<int> &inliers,
                                          const std::vector<double> &model_coefficients )
{
    // Get the plane normal
    // Calculate the 2-norm: norm (x) = sqrt (sum (abs (v)^2))
    // double n_norm = sqrt (model_coefficients.at (0) * model_coefficients.at (0) +
    //                      model_coefficients.at (1) * model_coefficients.at (1) +
    //                      model_coefficients.at (2) * model_coefficients.at (2));
    // model_coefficients.at (0) /= n_norm;
    // model_coefficients.at (1) /= n_norm;
    // model_coefficients.at (2) /= n_norm;
    // model_coefficients.at (3) /= n_norm;

    // Iterate through the 3d points and calculate the distances from them to the plane
    for( unsigned int i = 0; i < inliers.size( ); i++ )
    {
        // Calculate the distance from the point to the plane
        double distance_to_plane =
            model_coefficients.at( 0 ) * cloud_->points.at( inliers.at( i ) ).x +
            model_coefficients.at( 1 ) * cloud_->points.at( inliers.at( i ) ).y +
            model_coefficients.at( 2 ) * cloud_->points.at( inliers.at( i ) ).z +
            model_coefficients.at( 3 ) * 1;
        // Calculate the projection of the point on the plane
        cloud_->points.at( inliers.at( i ) ).x =
            cloud_->points.at( inliers.at( i ) ).x - distance_to_plane * model_coefficients.at( 0 );
        cloud_->points.at( inliers.at( i ) ).y =
            cloud_->points.at( inliers.at( i ) ).y - distance_to_plane * model_coefficients.at( 1 );
        cloud_->points.at( inliers.at( i ) ).z =
            cloud_->points.at( inliers.at( i ) ).z - distance_to_plane * model_coefficients.at( 2 );
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Check whether the given index samples can form a valid plane model, compute the model
 * coefficients from these samples and store them internally in model_coefficients_. The plane
 * coefficients are: a, b, c, d (ax+by+cz+d=0) \param samples the point indices found as possible
 * good candidates for creating a valid model
 */
bool SACModelPlane::computeModelCoefficients( const std::vector<int> &samples )
{
    model_coefficients_.resize( 4 );
    double Dx1, Dy1, Dz1, Dx2, Dy2, Dz2, Dy1Dy2;
    // Compute the segment values (in 3d) between XY
    Dx1 = cloud_->points.at( samples.at( 1 ) ).x - cloud_->points.at( samples.at( 0 ) ).x;
    Dy1 = cloud_->points.at( samples.at( 1 ) ).y - cloud_->points.at( samples.at( 0 ) ).y;
    Dz1 = cloud_->points.at( samples.at( 1 ) ).z - cloud_->points.at( samples.at( 0 ) ).z;

    // Compute the segment values (in 3d) between XZ
    Dx2 = cloud_->points.at( samples.at( 2 ) ).x - cloud_->points.at( samples.at( 0 ) ).x;
    Dy2 = cloud_->points.at( samples.at( 2 ) ).y - cloud_->points.at( samples.at( 0 ) ).y;
    Dz2 = cloud_->points.at( samples.at( 2 ) ).z - cloud_->points.at( samples.at( 0 ) ).z;

    Dy1Dy2 = Dy1 / Dy2;
    if( ( ( Dx1 / Dx2 ) == Dy1Dy2 ) && ( Dy1Dy2 == ( Dz1 / Dz2 ) ) ) // Check for collinearity
        return ( false );

    // Compute the plane coefficients from the 3 given points in a straightforward manner
    // calculate the plane normal n = (p2-p1) x (p3-p1) = cross (p2-p1, p3-p1)
    model_coefficients_[ 0 ] =
        ( cloud_->points.at( samples.at( 1 ) ).y - cloud_->points.at( samples.at( 0 ) ).y ) *
            ( cloud_->points.at( samples.at( 2 ) ).z - cloud_->points.at( samples.at( 0 ) ).z ) -
        ( cloud_->points.at( samples.at( 1 ) ).z - cloud_->points.at( samples.at( 0 ) ).z ) *
            ( cloud_->points.at( samples.at( 2 ) ).y - cloud_->points.at( samples.at( 0 ) ).y );

    model_coefficients_[ 1 ] =
        ( cloud_->points.at( samples.at( 1 ) ).z - cloud_->points.at( samples.at( 0 ) ).z ) *
            ( cloud_->points.at( samples.at( 2 ) ).x - cloud_->points.at( samples.at( 0 ) ).x ) -
        ( cloud_->points.at( samples.at( 1 ) ).x - cloud_->points.at( samples.at( 0 ) ).x ) *
            ( cloud_->points.at( samples.at( 2 ) ).z - cloud_->points.at( samples.at( 0 ) ).z );

    model_coefficients_[ 2 ] =
        ( cloud_->points.at( samples.at( 1 ) ).x - cloud_->points.at( samples.at( 0 ) ).x ) *
            ( cloud_->points.at( samples.at( 2 ) ).y - cloud_->points.at( samples.at( 0 ) ).y ) -
        ( cloud_->points.at( samples.at( 1 ) ).y - cloud_->points.at( samples.at( 0 ) ).y ) *
            ( cloud_->points.at( samples.at( 2 ) ).x - cloud_->points.at( samples.at( 0 ) ).x );
    // calculate the 2-norm: norm (x) = sqrt (sum (abs (v)^2))
    // nx ny nz (aka: ax + by + cz ...
    double n_norm = sqrt( model_coefficients_[ 0 ] * model_coefficients_[ 0 ] +
                          model_coefficients_[ 1 ] * model_coefficients_[ 1 ] +
                          model_coefficients_[ 2 ] * model_coefficients_[ 2 ] );
    model_coefficients_[ 0 ] /= n_norm;
    model_coefficients_[ 1 ] /= n_norm;
    model_coefficients_[ 2 ] /= n_norm;

    // ... + d = 0
    model_coefficients_[ 3 ] =
        -1 * ( model_coefficients_[ 0 ] * cloud_->points.at( samples.at( 0 ) ).x +
               model_coefficients_[ 1 ] * cloud_->points.at( samples.at( 0 ) ).y +
               model_coefficients_[ 2 ] * cloud_->points.at( samples.at( 0 ) ).z );

    return ( true );
}
#if 0
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /** \brief Recompute the plane coefficients using the given inlier set and return them to the user.
    * @note: these are the coefficients of the plane model after refinement (eg. after SVD)
    * \param inliers the data inliers found as supporting the model
    * \param refit_coefficients the resultant recomputed coefficients after non-linear optimization
    */
  void
    SACModelPlane::refitModel (const std::vector<int> &inliers, std::vector<double> &refit_coefficients)
  {
    if (inliers.size () == 0)
    {
      std::cout << "[SACModelPlane::RefitModel] Cannot re-fit 0 inliers!" << std::endl;
      refit_coefficients = model_coefficients_;
      return;
    }

    Eigen::Vector4d plane_coefficients;
    double curvature;

    // Use Least-Squares to fit the plane through all the given sample points and find out its coefficients
    cloud_geometry::nearest::computePointNormal (*cloud_, inliers, plane_coefficients, curvature);

    refit_coefficients.resize (4);
    for (int d = 0; d < 4; d++)
      refit_coefficients[d] = plane_coefficients (d);
  }
#endif
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Verify whether a subset of indices verifies the internal plane model coefficients.
 * \param indices the data indices that need to be tested against the plane model
 * \param threshold a maximum admissible distance threshold for determining the inliers from the
 * outliers
 */
bool SACModelPlane::doSamplesVerifyModel( const std::set<int> &indices, double threshold )
{
    for( std::set<int>::iterator it = indices.begin( ); it != indices.end( ); ++it )
        if( fabs( model_coefficients_.at( 0 ) * cloud_->points.at( *it ).x +
                  model_coefficients_.at( 1 ) * cloud_->points.at( *it ).y +
                  model_coefficients_.at( 2 ) * cloud_->points.at( *it ).z +
                  model_coefficients_.at( 3 ) ) > threshold )
            return ( false );

    return ( true );
}
} // namespace sample_consensus
