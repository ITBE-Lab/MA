##
# @package MABS.sweepAllReturnBest
# @brief Implements @ref MABS.sweepAllReturnBest.SweepAllReturnBest "SweepAllReturnBest".
# @file sweepAllReturnBest.py
# @brief Implements @ref MABS.sweepAllReturnBest.SweepAllReturnBest "SweepAllReturnBest".
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
    # Reimplemented from MABS.aligner.Module.get_input_type.
    def get_input_type(self):
        return [ContainerType.seedsVector, ContainerType.query, ContainerType.ref_seq]

    ##
    # @brief returns the @ref ContainerType "container type" seeds.
    # @details
    # Reimplemented from MABS.aligner.Module.get_output_type.
    def get_output_type(self):
        return ContainerType.seeds

    ##
    # @brief Execute LineSweep for all given seeds.
    # @details
    # Reimplemented from MABS.aligner.Module.execute.
    def execute(self, input):
        best_strips = SeedsVector()
        strips, query, ref_seq = input
        if len(strips) == 0:
            return Seeds()
        for strip in strips:
            app = self.__linesweep.execute((query, ref_seq, strip))
            best_strips.append(app)
        best = 0
        if len(best_strips) == 0:
            return Seeds()
        for index, strip in enumerate(best_strips):
            if strip.get_score() > best_strips[best].get_score():
                best = index
        #it's crucial that we return a new copy of the output since the input may be deleted
        return Seeds(best_strips[best])

