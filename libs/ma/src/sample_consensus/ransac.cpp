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
 * $Id: ransac.cpp 16379 2009-05-29 19:20:46Z hsujohnhsu $
 *
 */

/** \author Radu Bogdan Rusu */
#include <algorithm>
#include <cmath>
#include <limits>
#include <ma/sample_consensus/ransac.h>
#include <math.h>

namespace sample_consensus
{
////////////////////////////////////////////////////////////////////////////////
/** \brief RANSAC (RAndom SAmple Consensus) main constructor
 * \param model a Sample Consensus model
 * \param threshold distance to model threshold
 */
RANSAC::RANSAC( SACModel* model, double threshold ) : SAC( model )
{
    this->threshold_ = threshold;
    // Desired probability of choosing at least one sample free from outliers
    this->probability_ = 0.99;
    // Maximum number of trials before we give up.
    this->max_iterations_ = 100;

    this->iterations_ = 0;
}

////////////////////////////////////////////////////////////////////////////////
/** \brief RANSAC (RAndom SAmple Consensus) main constructor
 * \param model a Sample Consensus model
 */
RANSAC::RANSAC( SACModel* model ) : SAC( model )
{}

////////////////////////////////////////////////////////////////////////////////
/** \brief Compute the actual model and find the inliers
 * \param debug enable/disable on-screen debug information
 */
bool RANSAC::computeModel( int debug )
{
    iterations_ = 0;
    int n_best_inliers_count = -INT_MAX;
    double k = 1.0;

    std::vector<int> best_model;
    std::vector<int> best_inliers, inliers;
    std::vector<int> selection;

    int n_inliers_count = 0;

    // Iterate
    while( iterations_ < k )
    {
        // Get X samples which satisfy the model criteria
        sac_model_->getSamples( iterations_, selection );

        if( selection.size( ) == 0 )
            break;

        // Search for inliers in the point cloud for the current plane model M
        sac_model_->computeModelCoefficients( selection );
        /*
         * model coefficients are:
         * 0: x point 1
         * 1: y point 1
         * 2: z point 1 (unused)
         * 3: x point 2
         * 4: y point 2
         * 5: z point 2 (unused)
         */

        // added by markus
        // angle acceptable?
        double dHorizontalDist = sac_model_->getModelCoefficients( )[ 0 ] - sac_model_->getModelCoefficients( )[ 3 ];
        double dVerticalDist = sac_model_->getModelCoefficients( )[ 1 ] - sac_model_->getModelCoefficients( )[ 4 ];
        if( dHorizontalDist <= 0 && dVerticalDist <= 0 )
        {
            dHorizontalDist *= -1;
            dVerticalDist *= -1;
        } // if
        double dAngle = -90;
        // if we are in the correct quadrant we compute an angle
        if( dHorizontalDist > 0 && dVerticalDist > 0 )
            dAngle = atan( dVerticalDist / dHorizontalDist ) * 180 / std::acos( -1 );
        if( dAngle >= 20 && dAngle <= 70 )
        {

            sac_model_->selectWithinDistance( sac_model_->getModelCoefficients( ), threshold_, inliers );
            n_inliers_count = (int)inliers.size( );

            // Better match ?
            if( n_inliers_count > n_best_inliers_count )
            {
                // std::cout  << "New best angle: " << dAngle << std::endl;
                n_best_inliers_count = n_inliers_count;
                best_inliers = inliers;
                // inliers.clear ();
                best_model = selection;

                // Compute the k parameter (k=log(z)/log(1-w^n))
                double w = (double)( (double)n_inliers_count / (double)sac_model_->getIndices( )->size( ) );
                double p_no_outliers = 1 - pow( w, (double)selection.size( ) );
                p_no_outliers = std::max( std::numeric_limits<double>::epsilon( ),
                                          p_no_outliers ); // Avoid division by -Inf
                p_no_outliers = std::min( 1 - std::numeric_limits<double>::epsilon( ),
                                          p_no_outliers ); // Avoid division by 0.
                k = log( 1 - probability_ ) / log( p_no_outliers );
            }
        }
        else
            continue;

        iterations_ += 1;
        if( debug > 1 )
            std::cerr << "[RANSAC::computeModel] Trial " << iterations_ << " out of " << ceil( k ) << ": "
                      << n_inliers_count << " inliers (best is: " << n_best_inliers_count << " so far)." << std::endl;
        if( iterations_ > max_iterations_ )
        {
            if( debug > 0 )
                std::cerr << "[RANSAC::computeModel] RANSAC reached the maximum number of trials." << std::endl;
            break;
        }
    }

    if( best_model.size( ) != 0 )
    {
        if( debug > 0 )
            std::cerr << "[RANSAC::computeModel] Model found: " << n_best_inliers_count << " inliers." << std::endl;
        sac_model_->setBestModel( best_model );
        sac_model_->setBestInliers( best_inliers );
        return ( true );
    }
    else if( debug > 0 )
        std::cerr << "[RANSAC::computeModel] Unable to find a solution!" << std::endl;
    return ( false );
}
} // namespace sample_consensus
