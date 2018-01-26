##
# @package MA
# @file setupaligner.py
# @brief Implements @ref MA.set_up_aligner "set_up_aligner".
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
# - @ref MA.aligner.SweepAllReturnBest "SweepAllReturnBest"
# - NeedlemanWunsch
# - @ref MA.alignmentPrinter.AlignmentPrinter "AlignmentPrinter"
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
        num_strips=5,
        max_sweep = None,
        min_seeds= 2,
        min_seed_length= .05,
        max_seeds=7.0,
        nmw_give_up=20000
        ):

    soc = StripOfConsideration(max_hits, min_seeds, num_strips, min_seed_length, max_seeds)

    max_sweep_n = 0
    if not max_sweep is None:
        max_sweep_n = max_sweep

    execall = ExecOnVec(chain, not max_sweep is None, max_sweep_n)

    nmw = NeedlemanWunsch(nmw_give_up)
    nmw_multiple = ExecOnVec(nmw, True, 0)
    mappingQual = MappingQuality()

    extractAll = ExtractAllSeeds(max_hits)


    query_pledges_ = []
    if isinstance(query_pledges, list) or isinstance(query_pledges, tuple):
        query_pledges_.extend(query_pledges)
    else:
        query_pledges_.append(query_pledges)

    return_pledges = [[], [], [], [], []]

    for query_pledge in query_pledges_:
        ret_pl_indx = 0
        segment_pledge = seg.promise_me(
            fm_index_pledge,
            query_pledge
        )
        return_pledges[ret_pl_indx].append(segment_pledge)
        ret_pl_indx += 1

        strips_pledge = soc.promise_me(
            segment_pledge,
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

        align_pledge = mappingQual.promise_me(query_pledge, alignments_pledge)
        return_pledges[ret_pl_indx].append(align_pledge)
        ret_pl_indx += 1

    if isinstance(query_pledges, list) or isinstance(query_pledges, tuple):
        return return_pledges
    else:
        return [item for sublist in return_pledges for item in sublist]

def piped(
        query,
        reference,
        fm_index,
        seg=BinarySeeding(False),
        max_hits=5,
        num_strips=5,
        max_sweep = None,
        min_seeds= 2,
        min_seed_length= .05,
        max_seeds=7.0,
        nmw_give_up=20000,
        threads=32
        ):
    #static pledges
    fm_index_pledge = Pledge(FMIndex())
    fm_index_pledge.set(fm_index)

    reference_pledge = Pledge(Pack())
    reference_pledge.set(reference)

    #modules
    soc = StripOfConsideration(max_hits, min_seeds, num_strips, min_seed_length, max_seeds)

    max_sweep_n = 0
    if not max_sweep is None:
        max_sweep_n = max_sweep
    execall = ExecOnVec(chain, not max_sweep is None, max_sweep_n)

    nmw = NeedlemanWunsch(nmw_give_up)
    nmw_multiple = ExecOnVec(nmw, True, 0)
    getBestOnly = Tail(Alignment())

    #setup the pipes
    pipes = []
    for th in range(threads):
        #the in pledge for the pipe
        query_pledge = Pledge(NucSeq())
        segment_pledge = seg.promise_me(
            fm_index_pledge,
            query_pledge
        )
        strips_pledge = soc.promise_me(
            segment_pledge,
            query_pledge,
            reference_pledge,
            fm_index_pledge
        )

        chains_pledge = execall.promise_me(
            strips_pledge,
        )


        alignments_pledge = nmw_multiple.promise_me(
            chains_pledge,
            query_pledge,
            reference_pledge
        )
        #the out pledge for the pipe
        align_pledge = getBestOnly.promise_me(alignments_pledge)

        pipe = Pipe([query_pledge], align_pledge)
        pipes.append(pipe)

    #setup the query vector pledge
    query_vec_pledge = []
    for query in queries:
        query_vec_pledge.append(Pledge(NucSeq()))
        query_vec_pledge[-1].set(query)

    pipe_outs = []
    for pipe in pipes:
        pipe_outs.append(pipe.promise_me(query_vec_pledge))
    
    Pledge.simultaneous_get(pipe_outs, threads)

    out = []
    for pipe_out in pipe_outs:
        out.extend(pipe_out)

    return out