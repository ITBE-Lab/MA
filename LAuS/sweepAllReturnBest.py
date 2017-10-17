##
# @package LAuS.sweepAllReturnBest
# @brief Implements @ref LAuS.sweepAllReturnBest.SweepAllReturnBest "SweepAllReturnBest".
# @file sweepAllReturnBest.py
# @brief Implements @ref LAuS.sweepAllReturnBest.SweepAllReturnBest "SweepAllReturnBest".
# @author Markus Schmidt

from .aligner import *

##
# @brief Execute LineSweep for all given seeds.
# @details 
# The module calls LineSweep for all Seeds in the SeedsVector.
# @ingroup module
#
class SweepAllReturnBest(Module):

    def __init__(self):
        self.__linesweep = LineSweep()

    def __del__(self):
        print("SweepAllReturnBest destroyed")

    ##
    # @brief returns the @ref ContainerType "container types" seedsVector, query, ref_seq.
    # @details
    # Reimplemented from LAuS.aligner.Module.get_input_type.
    def get_input_type(self):
        return [ContainerType.seedsVector, ContainerType.query, ContainerType.ref_seq]

    ##
    # @brief returns the @ref ContainerType "container type" seeds.
    # @details
    # Reimplemented from LAuS.aligner.Module.get_output_type.
    def get_output_type(self):
        return ContainerType.seeds

    ##
    # @brief Execute LineSweep for all given seeds.
    # @details
    # Reimplemented from LAuS.aligner.Module.execute.
    def execute(self, input):
        print("PYTHON GOT CALLED")
        best_strips = []
        if input is None:
            print("whaaa")
        print("a")
        strips, query, ref_seq = input
        if len(strips) == 0:
            return Seeds()
        print("b...:")
        #print(str(strips[0] is None))
        print("b")
        for strip in strips:
            print("b.1")
            app = self.__linesweep.execute((query, ref_seq, strip))
            print("b.2")
            best_strips.append(app)
        
        print("c")

        best = 0
        for index, strip in enumerate(best_strips):
            if strip.get_score() > best_strips[best].get_score():
                best = index

        print("d")
        print("PYTHON RETURNING")
        #it's crucial that we return a new copy of the output since the input may be deleted
        return Seeds(best_strips[best])

