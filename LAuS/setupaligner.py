##
# @package LAuS
# @file setupaligner.py
# @brief Implements @ref LAuS.set_up_aligner "set_up_aligner".
# @author Markus Schmidt

from .__init__ import *

##
# @brief Setup the @ref comp_graph_sec "computational graph" 
# with the standard aligner comfiguration.
# @details 
# Uses following Modules in this order:
# - The given segmentation; default is LongestNonEnclosedSegments
# - NlongestIntervalsAsAnchors
# - Bucketing
# - @ref LAuS.aligner.SweepAllReturnBest "SweepAllReturnBest"
# - NeedlemanWunsch
# - @ref LAuS.alignmentPrinter.AlignmentPrinter "AlignmentPrinter"
#
# If a list of query pledges is given a list of result pledges is returned.
#
#
# Returns a list with one element for each query_pledge given.
# Each list element is a tuple of the returned pledges of every step in the alignment.
#
# @returns a list of pledge tuples
# @ingroup module
#
def set_up_aligner(query_pledges, reference_pledge, 
        fm_index_pledge, seg=LongestNonEnclosedSegments(), 
        chain=LineSweep(), max_hits=500, num_anchors=10):

    anc = NlongestIntervalsAsAnchors(num_anchors)

    bucketing = Bucketing()
    bucketing.max_hits = max_hits

    execall = SweepAllReturnBest(chain)

    nmw = NeedlemanWunsch()

    printer = AlignmentPrinter()

    query_pledges_ = []
    if isinstance(query_pledges, list) or isinstance(query_pledges, tuple):
        query_pledges_.extend(query_pledges)
    else:
        query_pledges_.append(query_pledges)

    return_pledges = ([], [], [], [], [])
    for query_pledge in query_pledges_:
        segment_pledge = seg.promise_me((
            fm_index_pledge,
            query_pledge
        ))
        return_pledges[0].append(segment_pledge)

        anchors_pledge = anc.promise_me((segment_pledge,))
        return_pledges[1].append(anchors_pledge)


        strips_pledge = bucketing.promise_me((
            segment_pledge,
            anchors_pledge,
            query_pledge,
            reference_pledge,
            fm_index_pledge
        ))
        return_pledges[2].append(strips_pledge)

        best_pledge = execall.promise_me((
            strips_pledge,
        ))
        return_pledges[3].append(best_pledge)

        align_pledge = nmw.promise_me((
            best_pledge,
            query_pledge,
            reference_pledge
        ))
        return_pledges[4].append(align_pledge)

    if isinstance(query_pledges, list) or isinstance(query_pledges, tuple):
        return return_pledges
    else:
        return return_pledges[0]