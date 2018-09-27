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
 * $Id: sac_model.h 21050 2009-08-07 21:24:30Z jfaustwg $
 *
 */

/** \author Radu Bogdan Rusu */

#ifndef _SAMPLE_CONSENSUS_SACMODEL_H_
#define _SAMPLE_CONSENSUS_SACMODEL_H_

//#include <geometry_msgs/Point32.h>     // ROS float point type
//#include <sensor_msgs/PointCloud.h>  // ROS point cloud type

#include <climits>
#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace geometry_msgs
{
struct Point32
{
    /*
     * This contains the position of a point in free space(with 32 bits of precision).
     * It is recommeded to use Point wherever possible instead of Point32.
     *
     * This recommendation is to promote interoperability.
     *
     * This message is designed to take up less space when sending
     * lots of points at once, as in the case of a PointCloud.
     */
    double x, y, z;
}; // struct
} // namespace geometry_msgs

namespace cloud_geometry
{
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Compute the cross product between two points (vectors).
 * \param p1 the first point/vector
 * \param p2 the second point/vector
 */
inline geometry_msgs::Point32 cross( const geometry_msgs::Point32 &p1,
                                     const geometry_msgs::Point32 &p2 )
{
    geometry_msgs::Point32 r;
    r.x = p1.y * p2.z - p1.z * p2.y;
    r.y = p1.z * p2.x - p1.x * p2.z;
    r.z = p1.x * p2.y - p1.y * p2.x;
    return ( r );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Compute the cross product between two points (vectors).
 * \param p1 the first point/vector
 * \param p2 the second point/vector
 */
inline geometry_msgs::Point32 cross( const std::vector<double> &p1,
                                     const geometry_msgs::Point32 &p2 )
{
    geometry_msgs::Point32 r;
    r.x = p1.at( 1 ) * p2.z - p1.at( 2 ) * p2.y;
    r.y = p1.at( 2 ) * p2.x - p1.at( 0 ) * p2.z;
    r.z = p1.at( 0 ) * p2.y - p1.at( 1 ) * p2.x;
    return ( r );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Compute the cross product between two points (vectors).
 * \param p1 the first point/vector
 * \param p2 the second point/vector
 */
inline geometry_msgs::Point32 cross( const geometry_msgs::Point32 &p1,
                                     const std::vector<double> &p2 )
{
    geometry_msgs::Point32 r;
    r.x = p2.at( 1 ) * p1.z - p2.at( 2 ) * p1.y;
    r.y = p2.at( 2 ) * p1.x - p2.at( 0 ) * p1.z;
    r.z = p2.at( 0 ) * p1.y - p2.at( 1 ) * p1.x;
    return ( r );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Compute the cross product between two points (vectors).
 * \param p1 the first point/vector
 * \param p2 the second point/vector
 * \param the resultant third vector p3 = p1 x p2
 */
inline void cross( const std::vector<double> &p1, const std::vector<double> &p2,
                   std::vector<double> &p3 )
{
    if ( p3.size( ) != 3 )
        p3.resize( 3 );
    p3[ 0 ] = p1.at( 1 ) * p2.at( 2 ) - p1.at( 2 ) * p2.at( 1 );
    p3[ 1 ] = p1.at( 2 ) * p2.at( 0 ) - p1.at( 0 ) * p2.at( 2 );
    p3[ 2 ] = p1.at( 0 ) * p2.at( 1 ) - p1.at( 1 ) * p2.at( 0 );
}
} // namespace cloud_geometry

namespace sensor_msgs
{
struct ChannelFloat32
{
    /*
     * This message is used by the PointCloud message to hold optional data
     * associated with each point in the cloud. The length of the values
     * array should be the same as the length of the points array in the
     * PointCloud, and each value should be associated with the corresponding
     * point.

     * Channel names in existing practice include:
     *   "u", "v" - row and column (respectively) in the left stereo image.
     *              This is opposite to usual conventions but remains for
     *              historical reasons. The newer PointCloud2 message has no
     *              such problem.
     *   "rgb" - For point clouds produced by color stereo cameras. uint8
     *           (R,G,B) values packed into the least significant 24 bits,
     *           in order.
     *   "intensity" - laser or pixel intensity.
     *   "distance"

     * The channel name should give semantics of the channel (e.g.
     * "intensity" instead of "value").
     */
    std::string name;
    /*
     * The values array should be 1-1 with the elements of the associated
     * PointCloud.
     */
    std::vector<double> values;
};

struct PointCloud
{
    /*
     * Array of 3d points. Each Point32 should be interpreted as a 3d point
     * in the frame given in the header.
     */
    std::vector<geometry_msgs::Point32> points;
    /*
     * Each channel should have the same number of elements as points array,
     * and the data in each channel should correspond 1:1 with each point.
     * Channel names in common practice are listed in ChannelFloat32.msg.
     */
    std::vector<ChannelFloat32> channels;

    void set_channels_size( unsigned int x )
    {
        channels.resize( x );
    }
    unsigned int get_channels_size( )
    {
        return channels.size( );
    }
}; // class
} // namespace sensor_msgs


namespace sample_consensus
{
class SACModel
{
   public:
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Constructor for base SACModel. */
    SACModel( ) : cloud_( NULL )
    {
    }
    SACModel( sensor_msgs::PointCloud cloud ) : cloud_( &cloud )
    {
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Destructor for base SACModel. */
    virtual ~SACModel( )
    {
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Get a set of random data samples and return them as point indices. Pure virtual.
     * \param iterations the internal number of iterations used by SAC methods
     * \param samples the resultant model samples
     */
    virtual void getSamples( int &iterations, std::vector<int> &samples ) = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Test whether the given model coefficients are valid given the input point cloud data.
     * Pure virtual. \param model_coefficients the model coefficients that need to be tested
     */
    virtual bool testModelCoefficients( const std::vector<double> &model_coefficients ) = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Check whether the given index samples can form a valid model, compute the model
     * coefficients from these samples and store them internally in model_coefficients_. Pure
     * virtual. \param samples the point indices found as possible good candidates for creating a
     * valid model
     */
    virtual bool computeModelCoefficients( const std::vector<int> &samples ) = 0;

#if 0
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** \brief Recompute the model coefficients using the given inlier set and return them to the user. Pure virtual.
        * @note: these are the coefficients of the model after refinement (eg. after a least-squares optimization)
        * \param inliers the data inliers found as supporting the model
        * \param refit_coefficients the resultant recomputed coefficients after non-linear optimization
        */
      virtual void refitModel (const std::vector<int> &inliers, std::vector<double> &refit_coefficients) = 0;
#endif
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Compute all distances from the cloud data to a given model. Pure virtual.
     * \param model_coefficients the coefficients of a model that we need to compute distances to
     * \param distances the resultant estimated distances
     */
    virtual void getDistancesToModel( const std::vector<double> &model_coefficients,
                                      std::vector<double> &distances ) = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Select all the points which respect the given model coefficients as inliers. Pure
     * virtual. \param model_coefficients the coefficients of a model that we need to compute
     * distances to \param threshold a maximum admissible distance threshold for determining the
     * inliers from the outliers \param inliers the resultant model inliers
     * @note: To get the refined inliers of a model, use:
     *        ANNpoint refined_coeff = refitModel (...); selectWithinDistance (refined_coeff,
     * threshold);
     */
    virtual void selectWithinDistance( const std::vector<double> &model_coefficients,
                                       double threshold, std::vector<int> &inliers ) = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Create a new point cloud with inliers projected onto the model. Pure virtual.
     * \param inliers the data inliers that we want to project on the model
     * \param model_coefficients the coefficients of a model
     * \param projected_points the resultant projected points
     */
    virtual void projectPoints( const std::vector<int> &inliers,
                                const std::vector<double> &model_coefficients,
                                sensor_msgs::PointCloud &projected_points ) = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Project inliers (in place) onto the given model. Pure virtual.
     * \param inliers the data inliers that we want to project on the model
     * \param model_coefficients the coefficients of a model
     */
    virtual void projectPointsInPlace( const std::vector<int> &inliers,
                                       const std::vector<double> &model_coefficients ) = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Verify whether a subset of indices verifies the internal model coefficients. Pure
     * virtual. \param indices the data indices that need to be tested against the model \param
     * threshold a maximum admissible distance threshold for determining the inliers from the
     * outliers
     */
    virtual bool doSamplesVerifyModel( const std::set<int> &indices, double threshold ) = 0;


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Set the dataset
     * \param cloud the data set to be used */
    inline void setDataSet( sensor_msgs::PointCloud *cloud )
    {
        this->cloud_ = cloud;
        indices_.clear( );
        indices_.resize( cloud_->points.size( ) );
        for ( unsigned int i = 0; i < cloud_->points.size( ); i++ )
            indices_[ i ] = i;
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Set the dataset and indices
     * \param cloud the data set to be used
     * \param indices the point indices used */
    inline void setDataSet( sensor_msgs::PointCloud *cloud, std::vector<int> indices )
    {
        this->cloud_ = cloud;
        this->indices_ = indices;
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Set the indices
     * \param indices the point indices used */
    void setDataIndices( std::vector<int> indices )
    {
        this->indices_ = indices;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Remove the inliers found from the initial set of given point indices. */
    virtual int removeInliers( );

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Return an unique id for each type of model employed. */
    virtual int getModelType( ) = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Set the best set of inliers. Used by SAC methods. Do not call this except if you know
     * what you're doing. \param best_inliers the set of inliers for the best model */
    void setBestInliers( const std::vector<int> &best_inliers )
    {
        this->best_inliers_ = best_inliers;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Return the best set of inliers found so far for this model. */
    std::vector<int> getBestInliers( )
    {
        return ( this->best_inliers_ );
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Set the best model. Used by SAC methods. Do not call this except if you know what
     * you're doing. \param best_model the best model found so far */
    void setBestModel( std::vector<int> best_model )
    {
        this->best_model_ = best_model;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Return the best model found so far. */
    std::vector<int> getBestModel( )
    {
        return ( this->best_model_ );
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** \brief Return the model coefficients of the best model found so far. */
    std::vector<double> getModelCoefficients( )
    {
        return ( this->model_coefficients_ );
    }

    /** \brief Return a pointer to the point cloud data. */
    sensor_msgs::PointCloud *getCloud( )
    {
        return ( this->cloud_ );
    }

    /** \brief Return a pointer to the point cloud data indices. */
    std::vector<int> *getIndices( )
    {
        return ( &this->indices_ );
    }

   protected:
    /** \brief Holds a pointer to the point cloud data array, since we don't want to copy the whole
     * thing here */
    sensor_msgs::PointCloud *cloud_;

    /** \brief The list of internal point indices used */
    std::vector<int> indices_;

    /** \brief The coefficients of our model computed directly from the best samples found */
    std::vector<double> model_coefficients_;

    /** \brief The model found after the last computeModel () as pointcloud indices */
    std::vector<int> best_model_;
    /** \brief The indices of the points that were chosen as inliers after the last computeModel ()
     * call */
    std::vector<int> best_inliers_;
};
} // namespace sample_consensus

#endif
