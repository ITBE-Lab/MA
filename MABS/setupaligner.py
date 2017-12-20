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
# - The given segmentation; default is BinarySeeding
# - GetAnchors
# - StripOfConsideration
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
        seg=BinarySeeding(False),
        chain=LinearLineSweep(),
        max_hits=5,
        num_anchors=5,
        strips_of_consideration=True,
        re_seed = None,
        max_sweep = None,
        strip_size = 1000
        ):

    anc = GetAnchors(num_anchors, max_hits)

    soc = StripOfConsideration(strip_size, max_hits)

    max_sweep_n = 0
    if not max_sweep is None:
        max_sweep_n = max_sweep

    execall = ExecOnVec(chain, not max_sweep is None, max_sweep_n)

    nmw = NeedlemanWunsch()
    nmw_multiple = ExecOnVec(nmw, True, 0)
    getBestOnly = Tail(Alignment())

    extractAll = ExtractAllSeeds(max_hits)

    reseed = ReSeed()
    if not re_seed is None:
        reseed.min_split_len = re_seed


    query_pledges_ = []
    if isinstance(query_pledges, list) or isinstance(query_pledges, tuple):
        query_pledges_.extend(query_pledges)
    else:
        query_pledges_.append(query_pledges)

    return_pledges = [[], [], [], []]
    if not re_seed is None:
        return_pledges.append([])
    if strips_of_consideration:
        return_pledges.append([])
        return_pledges.append([])

    for query_pledge in query_pledges_:
        ret_pl_indx = 0
        segment_pledge = seg.promise_me(
            fm_index_pledge,
            query_pledge
        )
        return_pledges[ret_pl_indx].append(segment_pledge)
        ret_pl_indx += 1

        if not re_seed is None:
            segment_pledge = reseed.promise_me(
                fm_index_pledge,
                segment_pledge,
                query_pledge
            )
            return_pledges[ret_pl_indx].append(segment_pledge)
            ret_pl_indx += 1


        if strips_of_consideration:
            anchors_pledge = anc.promise_me(segment_pledge,reference_pledge,fm_index_pledge)
            return_pledges[ret_pl_indx].append(anchors_pledge)
            ret_pl_indx += 1

            strips_pledge = soc.promise_me(
                segment_pledge,
                anchors_pledge,
                query_pledge,
                reference_pledge,
                fm_index_pledge
            )
            return_pledges[ret_pl_indx].append(strips_pledge)
            ret_pl_indx += 1

            chains_pledge = execall.promise_me(
                strips_pledge,
            )

            return_pledges[ret_pl_indx].append(chains_pledge)
            ret_pl_indx += 1

            alignments_pledge = nmw_multiple.promise_me(
                chains_pledge,
                query_pledge,
                reference_pledge
            )
            return_pledges[ret_pl_indx].append(alignments_pledge)
            ret_pl_indx += 1

            align_pledge = getBestOnly.promise_me(alignments_pledge)
            return_pledges[ret_pl_indx].append(align_pledge)
            ret_pl_indx += 1

        else:
            strip_pledge = extractAll.promise_me(
                segment_pledge, fm_index_pledge
            )
            return_pledges[ret_pl_indx].append(strip_pledge)
            ret_pl_indx += 1

            
            best_pledge = chain.promise_me(strip_pledge)

            return_pledges[ret_pl_indx].append(best_pledge)
            ret_pl_indx += 1

            align_pledge = nmw.promise_me(
                best_pledge,
                query_pledge,
                reference_pledge
            )
            return_pledges[ret_pl_indx].append(align_pledge)
            ret_pl_indx += 1

    if isinstance(query_pledges, list) or isinstance(query_pledges, tuple):
        return return_pledges
    else:
        return [item for sublist in return_pledges for item in sublist]