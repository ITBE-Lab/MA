##
# @package MA
# @file setupaligner.py
# @brief Implements @ref MA.set_up_aligner "set_up_aligner".
# @author Markus Schmidt

from .aligner import *

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
    ##
    # @brief sets up a new computational graph for a aligner
    #
    def __init__(
            self,
            max_hits=100,
            num_strips=5,
            complete_seeds = False,
            threads = 32,
            local = False
        ):
        self.query_vec_pledge = Pledge(ContainerVector(NucSeq()))
        self.reference_pledge = Pledge(Pack())
        self.fm_index_pledge = Pledge(FMIndex())
        self.__pledges = [ [], [], [], [] ]

        splitter = Splitter(self.query_vec_pledge)
        lock = Lock(NucSeq())
        seeding = BinarySeeding(complete_seeds)
        soc = StripOfConsideration(max_hits, num_strips)
        couple = ExecOnVec(LinearLineSweep())
        optimal = ExecOnVec(NeedlemanWunsch(local))
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
            self.__pledges[0].append(seeding_pledge)

            strips_pledge = soc.promise_me(
                seeding_pledge,
                query_pledge,
                self.reference_pledge,
                self.fm_index_pledge
            )
            self.__pledges[1].append(strips_pledge)

            couple_pledge = couple.promise_me(strips_pledge)
            self.__pledges[2].append(couple_pledge)

            alignments_pledge = optimal.promise_me(
                couple_pledge,
                query_pledge,
                self.reference_pledge
            )
            self.__pledges[3].append(alignments_pledge)

            align_pledge = mappingQual.promise_me(query_pledge, alignments_pledge)

            collector_pledge = self.collector.promise_me(align_pledge)

            unlock_pledge = unlock.promise_me(collector_pledge)

            self.return_pledges.append(unlock_pledge)

    ##
    # @brief sets the reference
    #
    def setRef(self, pack):
        self.reference_pledge.set(pack)

    ##
    # @brief sets the queries
    #
    def setQueries(self, queries):
        vec = ContainerVector(NucSeq())
        #@fixme this is due to a bug in the vec initialization...
        del vec[:]
        vec.extend(queries)
        self.query_vec_pledge.set(vec)

    ##
    # @brief sets the reference index
    #
    def setInd(self, index):
        self.fm_index_pledge.set(index)

    ##
    # @brief trigger the alignment optionally sets the queries
    #
    def align(self, queries = None):
        #reset runtimes
        for row in self.__pledges:
            for ele in row:
                ele.exec_time = 0
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

    ##
    # @brief returns the runtimes of the important stages
    #
    def get_runtimes(self):
        ret_list = {
            'seeding': 0,
            'seed processing: filtering': 0,
            'seed processing: coupling': 0,
            'optimal alignment': 0
        }
        for ele in self.__pledges[0]:
            ret_list['seeding'] += ele.exec_time
        for ele in self.__pledges[1]:
            ret_list['seed processing: filtering'] += ele.exec_time
        for ele in self.__pledges[2]:
            ret_list['seed processing: coupling'] += ele.exec_time
        for ele in self.__pledges[3]:
            ret_list['optimal alignment'] += ele.exec_time
        return ret_list



