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

import libLAuS
#TODO: remove this once unnecessary
from libLAuS import *
import traceback

##
# @brief the Python implementation of @ref Module "module".
# @details 
# The module overrides the @ref promiseMe "promise_me" function.
# Thus telling the C++ code that any module inheriting from this class is written in python.
# @see the C++ implementation of @ref Module "module".
# @ingroup module 
#
class Module(libLAuS.Module):

    ##
    # @brief The expected input types.
    # @details
    # Reimplemented from @ref Module::getInputType.
    def get_input_type(self):
        return [ContainerType.nothing]

    ##
    # @brief The expected output type.
    # @details
    # Reimplemented from @ref Module::getOutputType.
    def get_output_type(self):
        return ContainerType.nothing

    ##
    # @brief Execute the implemented algorithm.
    # @details
    # Reimplemented from @ref Module::saveExecute.
    def execute(self, input):
        self.__store_result = Module.execute(input)
        return self.__store_result

    ##
    # @brief call the execute function with a try catch block
    # @details
    # Reimplemented from @ref Module::saveExecute.
    def save_execute(self, input):
        try:
            return self.execute(input)
        except Exception as e:
            traceback.print_exc()
            return None

    ##
    # @brief Make this module promise to execute it's function on the provided data.
    # @details
    # Reimplemented from @ref Module::promiseMe.
    def promise_me(self, input):
        return Pledge.make_pledge(self, self.get_output_type(), input)

##
# @brief python wrapper for BinarySeeding
class BinarySeeding(libLAuS.BinarySeeding):
    def execute(self, *args):
        return super(BinarySeeding, self).execute(list(args))

    def promise_me(self, *args):
        return super(BinarySeeding, self).promise_me(list(args))

##
# @brief python wrapper for StripOfConsideration
class StripOfConsideration(libLAuS.StripOfConsideration):
    def execute(self, *args):
        return super(StripOfConsideration, self).execute(list(args))

    def promise_me(self, *args):
        return super(StripOfConsideration, self).promise_me(list(args))

##
# @brief python wrapper for LinearLineSweep
class LinearLineSweep(libLAuS.LinearLineSweep):
    def execute(self, *args):
        return super(LinearLineSweep, self).execute(list(args))

    def promise_me(self, *args):
        return super(LinearLineSweep, self).promise_me(list(args))

##
# @brief python wrapper for LinearLineSweep
class Chaining(libLAuS.Chaining):
    def execute(self, *args):
        return super(Chaining, self).execute(list(args))

    def promise_me(self, *args):
        return super(Chaining, self).promise_me(list(args))

##
# @brief python wrapper for LinearLineSweep
class NeedlemanWunsch(libLAuS.NeedlemanWunsch):
    def execute(self, *args):
        return super(NeedlemanWunsch, self).execute(list(args))

    def promise_me(self, *args):
        return super(NeedlemanWunsch, self).promise_me(list(args))

##
# @brief python wrapper for LinearLineSweep
class GetAnchors(libLAuS.GetAnchors):
    def execute(self, *args):
        return super(GetAnchors, self).execute(list(args))

    def promise_me(self, *args):
        return super(GetAnchors, self).promise_me(list(args))

"""
class ContainerType(CppContainerType):
    alignment = CppContainerType.alignment
    any = CppContainerType.any
    fM_index = CppContainerType.fM_index
    nothing = CppContainerType.nothing
    nucSeq = CppContainerType.nucSeq
    packedNucSeq = CppContainerType.packedNucSeq
    sa_interval = CppContainerType.sa_interval
    seed = CppContainerType.seed
    seeds = CppContainerType.seeds
    seedsVector = CppContainerType.seedsVector
    segment = CppContainerType.segment
    segmentList = CppContainerType.segmentList
    unknown = CppContainerType.unknown
    #from here the new python definitions start
    #whatever = len(CppContainerType.values)
"""
