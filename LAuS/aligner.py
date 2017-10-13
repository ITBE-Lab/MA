##
# @package LAuS
# @brief The Python part of the library is within this package
# @package LAuS.aligner
# @brief Implements @ref LAuS.aligner.Module "Module" and 
# @ref LAuS.aligner.SweepAllReturnBest "SweepAllReturnBest".
# @file aligner.py
# @brief Implements @ref LAuS.aligner.Module "Module" and 
# @ref LAuS.aligner.SweepAllReturnBest "SweepAllReturnBest".
# @author Markus Schmidt
#

from libLAuS import *

##
# @brief the Python implementation of @ref CppModule "module".
# @details 
# The module overrides the @ref promiseMe "promise_me" function.
# Thus telling the C++ code that any module inheriting from this class is written in python.
# @see the C++ implementation of @ref CppModule "module".
# @ingroup module 
#
class Module(CppModule):

    ##
    # @brief The expected input types.
    # @details
    # Reimplemented from @ref CppModule::getInputType.
    def get_input_type(self):
        return [ContainerType.nothing]

    ##
    # @brief The expected output type.
    # @details
    # Reimplemented from @ref CppModule::getOutputType.
    def get_output_type(self):
        return ContainerType.nothing

    ##
    # @brief Execute the implemented algorithm.
    # @details
    # Reimplemented from @ref CppModule::saveExecute.
    def execute(self, input):
        return Module.execute(input)

    ##
    # @brief Make this module promise to execute it's function on the provided data.
    # @details
    # Reimplemented from @ref CppModule::promiseMe.
    def promise_me(self, input):
        return Pledge(self, self.get_output_type(), input)

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

        return best_strip[best]

