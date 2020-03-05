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
 * $Id: sac.cpp 21050 2009-08-07 21:24:30Z jfaustwg $
 *
 */

/** \author Radu Bogdan Rusu */

#include <ma/sample_consensus/sac.h>

using namespace sample_consensus;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief return the point cloud representing a set of given indices.
 * \param indices a set of indices that represent the data that we're interested in */
sensor_msgs::PointCloud SAC::getPointCloud( std::vector<int> indices )
{
    sensor_msgs::PointCloud i_points;

    // Allocate enough space
    i_points.points.resize( indices.size( ) );
    i_points.set_channels_size( (unsigned int)sac_model_->getCloud( )->get_channels_size( ) );

    // Create the channels
    for( unsigned int d = 0; d < i_points.get_channels_size( ); d++ )
    {
        i_points.channels[ d ].name = sac_model_->getCloud( )->channels[ d ].name;
        i_points.channels[ d ].values.resize( indices.size( ) );
    }

    // Copy the data
    for( unsigned int i = 0; i < i_points.points.size( ); i++ )
    {
        i_points.points[ i ].x = sac_model_->getCloud( )->points[ indices.at( i ) ].x;
        i_points.points[ i ].y = sac_model_->getCloud( )->points[ indices.at( i ) ].y;
        i_points.points[ i ].z = sac_model_->getCloud( )->points[ indices.at( i ) ].z;
        for( unsigned int d = 0; d < i_points.get_channels_size( ); d++ )
            i_points.channels[ d ].values[ i ] = sac_model_->getCloud( )->channels[ d ].values[ indices.at( i ) ];
    }

    return ( i_points );
}
