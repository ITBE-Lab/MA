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
# - LongestNonEnclosedSegments
# - NlongestIntervalsAsAnchors
# - Bucketing
# - @ref LAuS.aligner.SweepAllReturnBest "SweepAllReturnBest"
# - NeedlemanWunsch
# - @ref LAuS.alignmentPrinter.AlignmentPrinter "AlignmentPrinter"
#
# If a list of query pledges is given a list of result pledges is returned.
#
# @returns the pledge of the alignment Printer.
# @ingroup module
#
def set_up_aligner(query_pledges, reference_pledge, fm_index_pledge):
    seg = LongestNonEnclosedSegments()

    anc = NlongestIntervalsAsAnchors(2)

    bucketing = Bucketing()

    sweep = SweepAllReturnBest()

    nmw = NeedlemanWunsch()

    printer = AlignmentPrinter()

    query_pledges_ = []
    if isinstance(query_pledges, list) or isinstance(query_pledges, tuple):
        query_pledges_.extend(query_pledges)
    else:
        query_pledges_.append(query_pledges)

    return_pledges = []
    for query_pledge in query_pledges_:
        segment_pledge = seg.promise_me((
            fm_index_pledge,
            query_pledge
        ))

        anchors_pledge = anc.promise_me((segment_pledge,))


        strips_pledge = bucketing.promise_me((
            segment_pledge,
            anchors_pledge,
            query_pledge,
            reference_pledge,
            fm_index_pledge
        ))

        best_pledge = sweep.promise_me((
            strips_pledge,
            query_pledge,
            reference_pledge
        ))

        align_pledge = nmw.promise_me((
            best_pledge,
            query_pledge,
            reference_pledge
        ))
        return_pledges.append(align_pledge)
        continue

        return_pledges.append(
            printer.promise_me((align_pledge, query_pledge, reference_pledge))
        )

    if isinstance(query_pledges, list) or isinstance(query_pledges, tuple):
        return return_pledges 
    else:
        return return_pledges[0]
