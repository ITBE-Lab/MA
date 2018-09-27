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
 * $Id: sac_model_line.cpp 21050 2009-08-07 21:24:30Z jfaustwg $
 *
 */

/** \author Radu Bogdan Rusu
 *
 * \todo Change the internal representation of the line model from 2 points to 1 point + direction.
 */

#include <sample_consensus/sac_model.h>
#include <sample_consensus/sac_model_line.h>
//#include <point_cloud_mapping/geometry/point.h>
//#include <point_cloud_mapping/geometry/nearest.h>

namespace sample_consensus
{
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Get 2 random points as data samples and return them as point indices.
 * \param iterations the internal number of iterations used by SAC methods
 * \param samples the resultant model samples
 * \note assumes unique points!
 */
void SACModelLine::getSamples( int &iterations, std::vector<int> &samples )
{
    samples.resize( 2 );
    double trand = indices_.size( ) / ( RAND_MAX + 1.0 );

    // Get a random number between 1 and max_indices
    int idx = ( int )( rand( ) * trand );
    // Get the index
    samples[ 0 ] = indices_.at( idx );

    // Get a second point which is different than the first
    int iter = 0;
    do
    {
        idx = ( int )( rand( ) * trand );
        samples[ 1 ] = indices_.at( idx );
        iter++;

        if ( iter > MAX_ITERATIONS_UNIQUE )
        {
            std::cout << "[SACModelLine::getSamples] WARNING: Could not select 2 unique points in "
                      << MAX_ITERATIONS_UNIQUE << " iterations!" << std::endl;
            break;
        }
        iterations++;
    } while ( samples[ 1 ] == samples[ 0 ] );
    iterations--;
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Select all the points which respect the given model coefficients as inliers.
 * \param model_coefficients the coefficients of a line model that we need to compute distances to
 * \param threshold a maximum admissible distance threshold for determining the inliers from the
 * outliers \param inliers the resultant model inliers \note: To get the refined inliers of a model,
 * use: ANNpoint refined_coeff = refitModel (...); selectWithinDistance (refined_coeff, threshold);
 */
void SACModelLine::selectWithinDistance( const std::vector<double> &model_coefficients,
                                         double threshold, std::vector<int> &inliers )
{
    double sqr_threshold = threshold * threshold;

    int nr_p = 0;
    inliers.resize( indices_.size( ) );

    // Obtain the line direction
    geometry_msgs::Point32 p3, p4;
    p3.x = model_coefficients.at( 3 ) - model_coefficients.at( 0 );
    p3.y = model_coefficients.at( 4 ) - model_coefficients.at( 1 );
    p3.z = model_coefficients.at( 5 ) - model_coefficients.at( 2 );

    // Iterate through the 3d points and calculate the distances from them to the plane
    for ( unsigned int i = 0; i < indices_.size( ); i++ )
    {
        // Calculate the distance from the point to the line
        // D = ||(P2-P1) x (P1-P0)|| / ||P2-P1|| = norm (cross (p2-p1, p2-p0)) / norm(p2-p1)
        // P1, P2 = line points, P0 = query point
        // P1P2 = <x2-x1, y2-y1, z2-z1> = <x3, y3, z3>
        // P1P0 = < x-x1,  y-y1, z-z1 > = <x4, y4, z4>
        // P1P2 x P1P0 = < y3*z4 - z3*y4, -(x3*z4 - x4*z3), x3*y4 - x4*y3 >
        //             = < (y2-y1)*(z-z1) - (z2-z1)*(y-y1), -[(x2-x1)*(z-z1) - (x-x1)*(z2-z1)],
        //             (x2-x1)*(y-y1) - (x-x1)*(y2-y1) >
        p4.x = model_coefficients.at( 3 ) - cloud_->points.at( indices_.at( i ) ).x;
        p4.y = model_coefficients.at( 4 ) - cloud_->points.at( indices_.at( i ) ).y;
        p4.z = model_coefficients.at( 5 ) - cloud_->points.at( indices_.at( i ) ).z;

        // P1P2 = sqrt (x3^2 + y3^2 + z3^2)
        // a = sqrt [(y3*z4 - z3*y4)^2 + (x3*z4 - x4*z3)^2 + (x3*y4 - x4*y3)^2]
        // double distance = SQR_NORM (cANN::cross (p4, p3)) / SQR_NORM (p3);
        geometry_msgs::Point32 c = cloud_geometry::cross( p4, p3 );
        double sqr_distance =
            ( c.x * c.x + c.y * c.y + c.z * c.z ) / ( p3.x * p3.x + p3.y * p3.y + p3.z * p3.z );

        if ( sqr_distance < sqr_threshold )
        {
            // Returns the indices of the points whose squared distances are smaller than the
            // threshold
            inliers[ nr_p ] = indices_[ i ];
            nr_p++;
        }
    }
    inliers.resize( nr_p );
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Compute all distances from the cloud data to a given line model.
 * \param model_coefficients the coefficients of a line model that we need to compute distances to
 * \param distances the resultant estimated distances
 */
void SACModelLine::getDistancesToModel( const std::vector<double> &model_coefficients,
                                        std::vector<double> &distances )
{
    distances.resize( indices_.size( ) );

    // Obtain the line direction
    geometry_msgs::Point32 p3, p4;
    p3.x = model_coefficients.at( 3 ) - model_coefficients.at( 0 );
    p3.y = model_coefficients.at( 4 ) - model_coefficients.at( 1 );
    p3.z = model_coefficients.at( 5 ) - model_coefficients.at( 2 );

    // Iterate through the 3d points and calculate the distances from them to the plane
    for ( unsigned int i = 0; i < indices_.size( ); i++ )
    {
        // Calculate the distance from the point to the line
        // D = ||(P2-P1) x (P1-P0)|| / ||P2-P1|| = norm (cross (p2-p1, p2-p0)) / norm(p2-p1)
        p4.x = model_coefficients.at( 3 ) - cloud_->points.at( indices_.at( i ) ).x;
        p4.y = model_coefficients.at( 4 ) - cloud_->points.at( indices_.at( i ) ).y;
        p4.z = model_coefficients.at( 5 ) - cloud_->points.at( indices_.at( i ) ).z;

        geometry_msgs::Point32 c = cloud_geometry::cross( p4, p3 );
        distances[ i ] =
            sqrt( c.x * c.x + c.y * c.y + c.z * c.z ) / ( p3.x * p3.x + p3.y * p3.y + p3.z * p3.z );
    }
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Create a new point cloud with inliers projected onto the line model.
 * \param inliers the data inliers that we want to project on the line model
 * \param model_coefficients the coefficients of a line model
 * \param projected_points the resultant projected points
 */
void SACModelLine::projectPoints( const std::vector<int> &inliers,
                                  const std::vector<double> &model_coefficients,
                                  sensor_msgs::PointCloud &projected_points )
{
    // Allocate enough space
    projected_points.points.resize( inliers.size( ) );
    projected_points.set_channels_size( cloud_->get_channels_size( ) );

    // Create the channels
    for ( unsigned int d = 0; d < projected_points.get_channels_size( ); d++ )
    {
        projected_points.channels[ d ].name = cloud_->channels[ d ].name;
        projected_points.channels[ d ].values.resize( inliers.size( ) );
    }

    // Compute the line direction (P2 - P1)
    geometry_msgs::Point32 p21;
    p21.x = model_coefficients.at( 3 ) - model_coefficients.at( 0 );
    p21.y = model_coefficients.at( 4 ) - model_coefficients.at( 1 );
    p21.z = model_coefficients.at( 5 ) - model_coefficients.at( 2 );

    // Iterate through the 3d points and calculate the distances from them to the line
    for ( unsigned int i = 0; i < inliers.size( ); i++ )
    {
        // double k = (DOT_PROD_3D (points[i], p21) - dotA_B) / dotB_B;
        double k = ( cloud_->points.at( inliers.at( i ) ).x * p21.x +
                     cloud_->points.at( inliers.at( i ) ).y * p21.y +
                     cloud_->points.at( inliers.at( i ) ).z * p21.z -
                     ( model_coefficients_[ 0 ] * p21.x + model_coefficients_[ 1 ] * p21.y +
                       model_coefficients_[ 2 ] * p21.z ) ) /
                   ( p21.x * p21.x + p21.y * p21.y + p21.z * p21.z );
        // Calculate the projection of the point on the line (pointProj = A + k * B)
        projected_points.points[ i ].x = model_coefficients_.at( 0 ) + k * p21.x;
        projected_points.points[ i ].y = model_coefficients_.at( 1 ) + k * p21.y;
        projected_points.points[ i ].z = model_coefficients_.at( 2 ) + k * p21.z;
        // Copy the other attributes
        for ( unsigned int d = 0; d < projected_points.get_channels_size( ); d++ )
            projected_points.channels[ d ].values[ i ] =
                cloud_->channels[ d ].values[ inliers.at( i ) ];
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Project inliers (in place) onto the given line model.
 * \param inliers the data inliers that we want to project on the line model
 * \param model_coefficients the coefficients of a line model
 */
void SACModelLine::projectPointsInPlace( const std::vector<int> &inliers,
                                         const std::vector<double> &model_coefficients )
{
    // Compute the line direction (P2 - P1)
    geometry_msgs::Point32 p21;
    p21.x = model_coefficients.at( 3 ) - model_coefficients.at( 0 );
    p21.y = model_coefficients.at( 4 ) - model_coefficients.at( 1 );
    p21.z = model_coefficients.at( 5 ) - model_coefficients.at( 2 );

    // Iterate through the 3d points and calculate the distances from them to the line
    for ( unsigned int i = 0; i < inliers.size( ); i++ )
    {
        // double k = (DOT_PROD_3D (points[i], p21) - dotA_B) / dotB_B;
        double k = ( cloud_->points.at( inliers.at( i ) ).x * p21.x +
                     cloud_->points.at( inliers.at( i ) ).y * p21.y +
                     cloud_->points.at( inliers.at( i ) ).z * p21.z -
                     ( model_coefficients_[ 0 ] * p21.x + model_coefficients_[ 1 ] * p21.y +
                       model_coefficients_[ 2 ] * p21.z ) ) /
                   ( p21.x * p21.x + p21.y * p21.y + p21.z * p21.z );
        // Calculate the projection of the point on the line (pointProj = A + k * B)
        cloud_->points.at( inliers.at( i ) ).x = model_coefficients_.at( 0 ) + k * p21.x;
        cloud_->points.at( inliers.at( i ) ).y = model_coefficients_.at( 1 ) + k * p21.y;
        cloud_->points.at( inliers.at( i ) ).z = model_coefficients_.at( 2 ) + k * p21.z;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Check whether the given index samples can form a valid line model, compute the model
 * coefficients from these samples and store them internally in model_coefficients_. The line
 * coefficients are represented by the points themselves. \param samples the point indices found as
 * possible good candidates for creating a valid model
 */
bool SACModelLine::computeModelCoefficients( const std::vector<int> &samples )
{
    model_coefficients_.resize( 6 );
    model_coefficients_[ 0 ] = cloud_->points.at( samples.at( 0 ) ).x;
    model_coefficients_[ 1 ] = cloud_->points.at( samples.at( 0 ) ).y;
    model_coefficients_[ 2 ] = cloud_->points.at( samples.at( 0 ) ).z;
    model_coefficients_[ 3 ] = cloud_->points.at( samples.at( 1 ) ).x;
    model_coefficients_[ 4 ] = cloud_->points.at( samples.at( 1 ) ).y;
    model_coefficients_[ 5 ] = cloud_->points.at( samples.at( 1 ) ).z;

    return ( true );
}

#if 0
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /** \brief Recompute the line coefficients using the given inlier set and return them to the user.
    * @note: these are the coefficients of the line model after refinement
    * \param inliers the data inliers found as supporting the model
    * \param refit_coefficients the resultant recomputed coefficients after non-linear optimization
    */
  void
    SACModelLine::refitModel (const std::vector<int> &inliers, std::vector<double> &refit_coefficients)
  {
    if (inliers.size () == 0)
    {
      std:cout << "[SACModelLine::RefitModel] Cannot re-fit 0 inliers!" << std::endl;
      refit_coefficients = model_coefficients_;
      return;
    }

    refit_coefficients.resize (6);

    // Compute the centroid of the samples
    geometry_msgs::Point32 centroid;
    // Compute the 3x3 covariance matrix
    Eigen::Matrix3d covariance_matrix;
    cloud_geometry::nearest::computeCovarianceMatrix (*cloud_, inliers, covariance_matrix, centroid);

    refit_coefficients[0] = centroid.x;
    refit_coefficients[1] = centroid.y;
    refit_coefficients[2] = centroid.z;

    // Extract the eigenvalues and eigenvectors
    //cloud_geometry::eigen_cov (covariance_matrix, eigen_values, eigen_vectors);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> ei_symm (covariance_matrix);
    Eigen::Vector3d eigen_values  = ei_symm.eigenvalues ();
    Eigen::Matrix3d eigen_vectors = ei_symm.eigenvectors ();

    refit_coefficients[3] = eigen_vectors (0, 2) + refit_coefficients[0];
    refit_coefficients[4] = eigen_vectors (1, 2) + refit_coefficients[1];
    refit_coefficients[5] = eigen_vectors (2, 2) + refit_coefficients[2];
  }
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Verify whether a subset of indices verifies the internal line model coefficients.
 * \param indices the data indices that need to be tested against the line model
 * \param threshold a maximum admissible distance threshold for determining the inliers from the
 * outliers
 */
bool SACModelLine::doSamplesVerifyModel( const std::set<int> &indices, double threshold )
{
    double sqr_threshold = threshold * threshold;
    for ( std::set<int>::iterator it = indices.begin( ); it != indices.end( ); ++it )
    {
        geometry_msgs::Point32 p3, p4;
        p3.x = model_coefficients_.at( 3 ) - model_coefficients_.at( 0 );
        p3.y = model_coefficients_.at( 4 ) - model_coefficients_.at( 1 );
        p3.z = model_coefficients_.at( 5 ) - model_coefficients_.at( 2 );

        p4.x = model_coefficients_.at( 3 ) - cloud_->points.at( *it ).x;
        p4.y = model_coefficients_.at( 4 ) - cloud_->points.at( *it ).y;
        p4.z = model_coefficients_.at( 5 ) - cloud_->points.at( *it ).z;

        geometry_msgs::Point32 c = cloud_geometry::cross( p4, p3 );
        double sqr_distance =
            ( c.x * c.x + c.y * c.y + c.z * c.z ) / ( p3.x * p3.x + p3.y * p3.y + p3.z * p3.z );

        if ( sqr_distance < sqr_threshold )
            return ( false );
    }
    return ( true );
}
} // namespace sample_consensus
