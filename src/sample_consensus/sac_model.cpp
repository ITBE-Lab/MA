/*
 * Copyright (c) 2008-2009 Radu Bogdan Rusu <rusu -=- cs.tum.edu>
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
 * $Id: sac_model.cpp 16379 2009-05-29 19:20:46Z hsujohnhsu $
 *
 */

/** \author Radu Bogdan Rusu */

#include <algorithm>
#include <iterator>
#include <sample_consensus/sac_model.h>

namespace sample_consensus
{
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Remove the inliers found from the initial set of given point indices. */
int SACModel::removeInliers( )
{
    // The point indices used for computing the current model are in indices_
    // What we need to do is subtract the Inliers
    std::vector<int> remaining_indices;

    // Sort the inliers and the point cloud indices
    std::sort( best_inliers_.begin( ), best_inliers_.end( ) );
    std::sort( indices_.begin( ), indices_.end( ) );

    set_difference( indices_.begin( ), indices_.end( ), best_inliers_.begin( ),
                    best_inliers_.end( ),
                    std::inserter( remaining_indices, remaining_indices.begin( ) ) );

    indices_ = remaining_indices;


    return indices_.size( );
}
} // namespace sample_consensus
