##
# @package MA.setupaligner
# @brief Implements @ref MA.set_up_aligner "set_up_aligner".
# @file setupaligner.py
# @brief Implements @ref MA.set_up_aligner "set_up_aligner".
# @author Markus Schmidt

from .aligner import *
from .postgresInterface import *

class ExecOnVecPy(Module):
    def __init__(self, module):
        self.module = module

    def get_input_type(self):
        module_input = self.module.get_input_type()
        module_input[0] = ContainerVector(module_input[0])
        return module_input

    def get_output_type(self):
        return ContainerVector(self.module.get_output_type())

    def execute(self, *input):
        ret = self.get_output_type()
        element_list = input[0]
        rem_inputs = input[1:]
        for element in element_list:
            # star does unpack the list into arguments (called "splat")
            ret.append(self.module.execute(element, *rem_inputs))
        return ret

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
    def __init__(self):
        self.query_vec_pledge = Pledge(ContainerVector(NucSeq()))
        self.ref_pledge = Pledge(Pack())
        self.fm_pledge = Pledge(FMIndex())

        splitter = Splitter(self.query_vec_pledge)
        lock = Lock(Nil())
        seeding = BinarySeeding()
        soc = StripOfConsideration()
        harmonization = LinearLineSweep()
        optimal = ExecOnVec(NeedlemanWunsch())
        mappingQual = MappingQuality()
        self.dbWriter = DbWriter()
        dbWriterLoop = ExecOnVecPy(self.dbWriter)

        non_existant_pledge = Pledge( Nil() )
        non_existant_pledge.set(Nil())

        self.query_pledge = lock.promise_me( splitter.promise_me( non_existant_pledge ) )
        unlock = UnLock(self.query_pledge)

        self.seed_pledge = seeding.promise_me( self.fm_pledge, self.query_pledge )
        self.soc_pledge = soc.promise_me(
                self.seed_pledge, self.query_pledge, self.ref_pledge, self.fm_pledge
            )
        self.harm_pledge = harmonization.promise_me(
                self.soc_pledge, self.query_pledge
            )
        self.nw_pledge = optimal.promise_me(
                self.harm_pledge, self.query_pledge, self.ref_pledge
            )
        self.map_pledge = mappingQual.promise_me(
            self.query_pledge, self.nw_pledge
        )
        self.writer_pledge = dbWriterLoop.promise_me(
            self.map_pledge, self.query_pledge,  self.ref_pledge
        )
        self.unlock_pledge = unlock.promise_me(self.writer_pledge)

    ##
    # @brief sets the reference
    #
    def setRef(self, pack):
        self.ref_pledge.set(pack)

    ##
    # @brief sets the queries
    #
    def setQueries(self, queries):
        vec = ContainerVector(NucSeq())
        vec.extend(queries)
        self.query_vec_pledge.set(vec)

    ##
    # @brief sets the reference index
    #
    def setInd(self, index):
        self.fm_pledge.set(index)

    ##
    # @brief trigger the alignment optionally sets the queries
    #
    ## def align(self, queries = None):
    ##     #reset runtimes
    ##     for row in self.__pledges:
    ##         for ele in row:
    ##             ele.exec_time = 0
    ##     if not queries is None:
    ##         self.setQueries(queries)
    ##     #the actual alignment
    ##     Pledge.simultaneous_get(self.return_pledges, True)
    ##     alignments = []
    ##     for alignment in self.collector.content:
    ##         alignments.append(alignment[0])
    ##     #remove the content
    ##     del self.collector.content[:]
    ##     return alignments
    def get_one(self):
        self.unlock_pledge.get()
        self.dbWriter.finalize()




