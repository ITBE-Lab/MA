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
def set_up_aligner(
        query_pledges, 
        reference_pledge, 
        fm_index_pledge, 
        seg=LongestNonEnclosedSegments(), 
        chain=LineSweep(), 
        max_hits=5, 
        num_anchors=10, 
        strips_of_consideration=True
        ):

    anc = NlongestIntervalsAsAnchors(num_anchors)

    bucketing = Bucketing()
    bucketing.max_hits = max_hits

    execall = SweepAllReturnBest(chain)

    nmw = NeedlemanWunsch()
    nmw_multiple = NmwMultiple()
    nmw_multiple.try_n_many = num_anchors
    getBestOnly = GetBestOnly()

    printer = AlignmentPrinter()

    extractAll = ExtractAllSeeds()
    extractAll.max_hits = max_hits

    query_pledges_ = []
    if isinstance(query_pledges, list) or isinstance(query_pledges, tuple):
        query_pledges_.extend(query_pledges)
    else:
        query_pledges_.append(query_pledges)

    return_pledges = [[], [], [], []]
    if strips_of_consideration:
        return_pledges.append([])

    for query_pledge in query_pledges_:
        ret_pl_indx = 0
        segment_pledge = seg.promise_me((
            fm_index_pledge,
            query_pledge
        ))
        return_pledges[ret_pl_indx].append(segment_pledge)
        ret_pl_indx += 1


        if strips_of_consideration:
            anchors_pledge = anc.promise_me((segment_pledge,))
            return_pledges[ret_pl_indx].append(anchors_pledge)
            ret_pl_indx += 1


            strips_pledge = bucketing.promise_me((
                segment_pledge,
                anchors_pledge,
                query_pledge,
                reference_pledge,
                fm_index_pledge
            ))
            return_pledges[ret_pl_indx].append(strips_pledge)
            ret_pl_indx += 1

            best_pledge = execall.promise_me((
                strips_pledge,
            ))

            return_pledges[ret_pl_indx].append(best_pledge)
            ret_pl_indx += 1

            align_pledge = nmw_multiple.promise_me((
                best_pledge,
                query_pledge,
                reference_pledge
            ))

            return_pledges[ret_pl_indx].append(align_pledge)
            ret_pl_indx += 1

        else:
            strip_pledge = extractAll.promise_me((
                segment_pledge, fm_index_pledge
            ))
            return_pledges[ret_pl_indx].append(strip_pledge)
            ret_pl_indx += 1

            
            best_pledge = chain.promise_me((
                strip_pledge,
            ))

            return_pledges[ret_pl_indx].append(best_pledge)
            ret_pl_indx += 1

            align_pledge = nmw.promise_me((
                best_pledge,
                query_pledge,
                reference_pledge
            ))
            return_pledges[ret_pl_indx].append(align_pledge)
            ret_pl_indx += 1

    if isinstance(query_pledges, list) or isinstance(query_pledges, tuple):
        return return_pledges
    else:
        return return_pledges[0]