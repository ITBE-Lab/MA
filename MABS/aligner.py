##
# @package MABS
# @brief The Python part of the library is within this package
# @package MABS.aligner
# @brief Implements @ref MABS.aligner.Module "Module" and 
# @ref MABS.aligner.SweepAllReturnBest "SweepAllReturnBest".
# @file aligner.py
# @brief Implements @ref MABS.aligner.Module "Module" and 
# @ref MABS.aligner.SweepAllReturnBest "SweepAllReturnBest".
# @author Markus Schmidt
#

import libMABS
import traceback

##
# @brief the Python implementation of @ref Module "module".
# @details 
# The module overrides the @ref promiseMe "promise_me" function.
# Thus telling the C++ code that any module inheriting from this class is written in python.
# @see the C++ implementation of @ref Module "module".
# @ingroup module 
#
class Module(libMABS.Module):

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
# @brief contains the final output of the aligner.
# @details
# Holds a sparse vector like representation of one alignment.
#
# @ingroup container
#
class Alignment(libMABS.Alignment):
    pass

##
# @brief Describes the type of match at one specific position of the alignment.
# @details
# - match: query and reference have the same nucleotide.
# - seed: query and reference have the same nucleotide (the match was found as part of a seed).
# - missmatch: query and reference have different nucleotides, 
#   but they are aligned to the same position nonetheless.
# - insertion: a nucleotide is present on the query that has no counterpart on the reference.
# - deletion: a nucleotide is present on the reference that has no counterpart on the query.
# @note libMABS::MatchType is an enum and therefore does not show up in the class hierarchy.
#
class MatchType(libMABS.MatchType):
    pass

##
# @brief Represents the pledge to deliver some container.
# @details
# Content may be provided by a @ref Module::promiseMe "module" or by calling @ref set.
#
# @ingroup container
#
class Pledge(libMABS.Pledge):
    pass

##
# @brief Contains a suffix array.
# @details
#
# @ingroup container
#
class FMIndex(libMABS.FMIndex):
    pass

##
# @brief And empty Container.
# @details
#
# @ingroup container
#
class Nil(libMABS.Nil):
    pass

##
# @brief The Container for a SegmentVector
# @details
# A doubly linked list holding Segments.
#
# @ingroup container
#
class SegmentVector(libMABS.SegmentVector):
    pass

##
# @brief The Container for a Seeds_
# @details
# is iterable
#
# @ingroup container
#
class Seeds(libMABS.Seeds):
    pass

##
# @brief A interval on the query that contains one or multiple @ref Seed "seeds".
# @details
#
# @ingroup container
#
class Segment(libMABS.Segment):
    pass

##
# @brief A packed version of a @ref NucSeq_ "NucSeq".
# @details
#
# @ingroup container
#
class Pack(libMABS.Pack):
    pass

##
# @brief A single seed.
# @details
#
# @ingroup container
#
class Seed(libMABS.Seed):
    pass

##
# @brief A SAInterval.
# @details
#
# @ingroup container
#
class SAInterval(libMABS.SAInterval):
    pass

##
# @brief A nucleotide sequence.
# @details
#
# @ingroup container
#
class NucSeq(libMABS.NucSeq):
    pass

##
# @brief The ContainerVector Module.
# @details
# Holds multiple containers.
#
# @ingroup container
#
class ContainerVector(libMABS.ContainerVector):
    def __init__(self):
        libMABS.ContainerVector.__init__(self)

    def __init__(self, *args):
        libMABS.ContainerVector.__init__(self, args[0])
        for arg in args:
            self.append(arg)

##
# @brief python wrapper for BinarySeeding
class BinarySeeding(libMABS.BinarySeeding):
    def execute(self, *args):
        return super(BinarySeeding, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(BinarySeeding, self).promise_me(ContainerVector(*args))

##
# @brief python wrapper for StripOfConsideration
class StripOfConsideration(libMABS.StripOfConsideration):
    def execute(self, *args):
        return super(StripOfConsideration, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(StripOfConsideration, self).promise_me(ContainerVector(*args))

##
# @brief python wrapper for LinearLineSweep
class LinearLineSweep(libMABS.LinearLineSweep):
    def execute(self, *args):
        return super(LinearLineSweep, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(LinearLineSweep, self).promise_me(ContainerVector(*args))

##
# @brief python wrapper for LinearLineSweep
class Chaining(libMABS.Chaining):
    def execute(self, *args):
        return super(Chaining, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(Chaining, self).promise_me(ContainerVector(*args))

##
# @brief python wrapper for LinearLineSweep
class NeedlemanWunsch(libMABS.NeedlemanWunsch):
    def execute(self, *args):
        return super(NeedlemanWunsch, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(NeedlemanWunsch, self).promise_me(ContainerVector(*args))

##
# @brief python wrapper for ExecOnVec
class ExecOnVec(libMABS.ExecOnVec):
    def execute(self, *args):
        return super(ExecOnVec, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(ExecOnVec, self).promise_me(ContainerVector(*args))

##
# @brief The Tail Module.
# @details
# returns the tail of a container vector
# @ingroup module
#
class Tail(libMABS.Tail):
    def execute(self, *args):
        return super(Tail, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(Tail, self).promise_me(ContainerVector(*args))

##
# @brief The ExtractAllSeeds Module.
# @details
# extracts all seeds from a segment list.
# @ingroup module
#
class ExtractAllSeeds(libMABS.ExtractAllSeeds):
    def execute(self, *args):
        return super(ExtractAllSeeds, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(ExtractAllSeeds, self).promise_me(ContainerVector(*args))

##
# @brief The ReSeed Module.
# @details
# extracts all seeds from a segment list.
# @ingroup module
#
class ReSeed(libMABS.ReSeed):
    def execute(self, *args):
        return super(ReSeed, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(ReSeed, self).promise_me(ContainerVector(*args))

##
# @brief The SMW Module.
# @ingroup module
#
class SMW(libMABS.SMW):
    def execute(self, *args):
        return super(SMW, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(SMW, self).promise_me(ContainerVector(*args))