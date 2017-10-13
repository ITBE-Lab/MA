##
# @package LAuS
# @file setupaligner.py
# @brief Implements @ref LAuS.set_up_aligner "set_up_aligner".
# @author Markus Schmidt

from .aligner import *

##
# @brief Setup the @ref comp_graph_sec "computational graph" 
# with the standard aligner comfiguration.
# @details 
# Uses following Modules in this order:
# - Segmentation
# - NlongestIntervalsAsAnchors
# - Bucketing
# - @ref LAuS.aligner.SweepAllReturnBest "SweepAllReturnBest"
# - NeedlemanWunsch
# - @ref LAuS.alignmentPrinter.AlignmentPrinter "AlignmentPrinter"
#
# @returns the pledge of the alignment Printer.
# @ingroup module
#
def set_up_aligner(query_pledge, reference_pledge, fm_index_pledge):
    seg = Segmentation()

    anc = NlongestIntervalsAsAnchors(2)

    bucketing = Bucketing()
    bucketing.strip_size = 50

    sweep = SweepAllReturnBest()

    nmw = NeedlemanWunsch()

    printer = AlignmentPrinter()


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

    return printer.promise_me((align_pledge, query_pledge, reference_pledge))