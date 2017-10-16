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
        strips, query, ref_seq = input
        best_strip = []
        for strip in strips:
            app = self.__linesweep.execute((query, ref_seq, strip))
            best_strip.append(app)

        best = 0
        for index, strip in enumerate(best_strip):
            if strip.get_score() > best_strip[best].get_score():
                best = index

        print("blub")
        return best_strip[best]

