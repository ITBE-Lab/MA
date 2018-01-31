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
class Aligner:
    def __init__(
                self,
                max_hits=5,
                num_strips=5,
                complete_seeds = False,
                threads = 32
            ):
        self.query_vec_pledge = Pledge(ContainerVector(NucSeq()))
        self.reference_pledge = Pledge(Pack())
        self.fm_index_pledge = Pledge(FMIndex())

        splitter = Splitter(self.query_vec_pledge)
        lock = Lock(NucSeq())
        seeding = BinarySeeding(complete_seeds)
        soc = StripOfConsideration(max_hits, num_strips)
        couple = ExecOnVec(LinearLineSweep())
        optimal = ExecOnVec(NeedlemanWunsch())
        mappingQual = MappingQuality()

        self.collector = Collector(NucSeq())
        self.return_pledges = []

        for _ in range(threads):
            ret_pl_indx = 0

            nil_pledge = Pledge(Nil())
            nil_pledge.set(Nil())
            query_pledge = lock.promise_me(
                splitter.promise_me(nil_pledge)
            )
            unlock = UnLock(query_pledge)

            seeding_pledge = seeding.promise_me(self.fm_index_pledge,query_pledge)

            strips_pledge = soc.promise_me(
                seeding_pledge,
                query_pledge,
                self.reference_pledge,
                self.fm_index_pledge
            )

            couple_pledge = couple.promise_me(strips_pledge)

            alignments_pledge = optimal.promise_me(
                couple_pledge,
                query_pledge,
                self.reference_pledge
            )

            align_pledge = mappingQual.promise_me(query_pledge, alignments_pledge)

            collector_pledge = self.collector.promise_me(align_pledge)

            unlock_pledge = unlock.promise_me(collector_pledge)

            self.return_pledges.append(unlock_pledge)

    def setRef(self, pack):
        self.reference_pledge.set(pack)

    def setQueries(self, queries):
        vec = ContainerVector(NucSeq())
        del vec[:]
        vec.extend(queries)
        self.query_vec_pledge.set(vec)

    def setRef(self, pack):
        self.reference_pledge.set(pack)

    def setInd(self, index):
        self.fm_index_pledge.set(index)

    def align(self, queries = None):
        if not queries is None:
            self.setQueries(queries)
        #the actual alignment
        Pledge.simultaneous_get(self.return_pledges, True)
        alignments = []
        for alignment in self.collector.content:
            alignments.append(alignment[0])
        #remove the content
        del self.collector.content[:]
        return alignments

