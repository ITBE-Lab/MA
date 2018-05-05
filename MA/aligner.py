##
# @package MA
# @brief The Python part of the library is within this package
# @package MA.aligner
# @brief Implements @ref MA.aligner.Module "Module" and 
# @ref MA.aligner.SweepAllReturnBest "SweepAllReturnBest".
# @file aligner.py
# @brief Implements @ref MA.aligner.Module "Module" and 
# @ref MA.aligner.SweepAllReturnBest "SweepAllReturnBest".
# @author Markus Schmidt
#

import libMA
import traceback

##
# @brief the Python implementation of @ref Module "module".
# @details 
# The module overrides the @ref promiseMe "promise_me" function.
# Thus telling the C++ code that any module inheriting from this class is written in python.
# @see the C++ implementation of @ref Module "module".
# @ingroup module 
#
class Module(libMA.Module):

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
    def execute(self, *args):
        self.__store_result = Module.execute(args)
        return self.__store_result

    ##
    # @brief call the execute function with a try catch block
    # @details
    # Reimplemented from @ref Module::saveExecute.
    def save_execute(self, args):
        try:
            #we get a tuple (args) from the cpp code that we then expand into multiple arguments so that execute can be called with an argument list
            return self.execute(*args)
        except Exception as e:
            traceback.print_exc()
            return None

    ##
    # @brief Make this module promise to execute it's function on the provided data.
    # @details
    # Reimplemented from @ref Module::promiseMe.
    def promise_me(self, *args):
        return Pledge.make_pledge(self, self.get_output_type(), args)


##
# @brief contains the final output of the aligner.
# @details
# Holds a sparse vector like representation of one alignment.
#
# @ingroup container
#
class Alignment(libMA.Alignment):
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
# @note libMA::MatchType is an enum and therefore does not show up in the class hierarchy.
#
class MatchType(libMA.MatchType):
    pass

##
# @brief Represents the pledge to deliver some container.
# @details
# Content may be provided by a @ref Module::promiseMe "module" or by calling @ref set.
#
# @ingroup container
#
class Pledge(libMA.Pledge):
    pass

##
# @brief Contains a suffix array.
# @details
#
# @ingroup container
#
class FMIndex(libMA.FMIndex):
    pass

##
# @brief And empty Container.
# @details
#
# @ingroup container
#
class Nil(libMA.Nil):
    pass

##
# @brief The Container for a SegmentVector
# @details
# A doubly linked list holding Segments.
#
# @ingroup container
#
class SegmentVector(libMA.SegmentVector):
    pass

##
# @brief The Container for a Seeds_
# @details
# is iterable
#
# @ingroup container
#
class Seeds(libMA.Seeds):
    pass

##
# @brief A interval on the query that contains one or multiple @ref Seed "seeds".
# @details
#
# @ingroup container
#
class Segment(libMA.Segment):
    pass

##
# @brief A packed version of a @ref NucSeq_ "NucSeq".
# @details
#
# @ingroup container
#
class Pack(libMA.Pack):
    pass

##
# @brief A single seed.
# @details
#
# @ingroup container
#
class Seed(libMA.Seed):
    pass

##
# @brief A SAInterval.
# @details
#
# @ingroup container
#
class SAInterval(libMA.SAInterval):
    pass

##
# @brief A nucleotide sequence.
# @details
#
# @ingroup container
#
class NucSeq(libMA.NucSeq):
    ##
    # @brief enable slicing a NucSeq
    def __getitem__(self, val):
        if isinstance(val, slice):
            ret = NucSeq()
            for index in range(val.start or 0, val.stop or len(self), val.step or 1):
                ret.append(libMA.NucSeq.__getitem__(self, index))
            return ret
        else:
            return libMA.NucSeq.__getitem__(self, val)

##
# @brief The ContainerVector Module.
# @details
# Holds multiple containers.
#
# @ingroup container
#
class ContainerVector(libMA.ContainerVector):
    def __init__(self, *args):
        if len(args) == 0:
            libMA.ContainerVector.__init__(self)
            return
        libMA.ContainerVector.__init__(self, args[0])
        for arg in args:
            self.append(arg)

##
# @brief python wrapper for BinarySeeding
class BinarySeeding(libMA.BinarySeeding):
    def execute(self, *args):
        return super(BinarySeeding, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(BinarySeeding, self).promise_me(ContainerVector(*args))

##
# @brief python wrapper for StripOfConsideration
class StripOfConsideration(libMA.StripOfConsideration):
    def execute(self, *args):
        return super(StripOfConsideration, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(StripOfConsideration, self).promise_me(ContainerVector(*args))


##
# @brief python wrapper for LinearLineSweep
class LinearLineSweep(libMA.LinearLineSweep):
    def execute(self, *args):
        return super(LinearLineSweep, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(LinearLineSweep, self).promise_me(ContainerVector(*args))

##
# @brief python wrapper for FileReader
class FileReader(libMA.FileReader):
    def execute(self):
        return super(FileReader, self).execute(ContainerVector(Nil()))

    def promise_me(self, *args):
        return super(FileReader, self).promise_me(ContainerVector(*args))

##
# @brief python wrapper for LinearLineSweep
class NeedlemanWunsch(libMA.NeedlemanWunsch):
    def execute(self, *args):
        return super(NeedlemanWunsch, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(NeedlemanWunsch, self).promise_me(ContainerVector(*args))

##
# @brief python wrapper for ExecOnVec
class ExecOnVec(libMA.ExecOnVec):
    def execute(self, *args):
        return super(ExecOnVec, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(ExecOnVec, self).promise_me(ContainerVector(*args))

##
# @brief python wrapper for MappingQuality
class MappingQuality(libMA.MappingQuality):
    def execute(self, *args):
        return super(MappingQuality, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(MappingQuality, self).promise_me(ContainerVector(*args))

##
# @brief The Tail Module.
# @details
# returns the tail of a container vector
# @ingroup module
#
class Tail(libMA.Tail):
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
class ExtractAllSeeds(libMA.ExtractAllSeeds):
    def execute(self, *args):
        return super(ExtractAllSeeds, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(ExtractAllSeeds, self).promise_me(ContainerVector(*args))

##
# @brief The Splitter Module.
# @ingroup module
#
class Splitter(libMA.Splitter):
    def execute(self, *args):
        return super(Splitter, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(Splitter, self).promise_me(ContainerVector(*args))

##
# @brief The Collector Module.
# @ingroup module
#
class Collector(libMA.Collector):
    def execute(self, *args):
        return super(Collector, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(Collector, self).promise_me(ContainerVector(*args))

##
# @brief The Lock Module.
# @ingroup module
#
class Lock(libMA.Lock):
    def execute(self, *args):
        return super(Lock, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(Lock, self).promise_me(ContainerVector(*args))

##
# @brief The UnLock Module.
# @ingroup module
#
class UnLock(libMA.UnLock):
    def execute(self, *args):
        return super(UnLock, self).execute(ContainerVector(*args))

    def promise_me(self, *args):
        return super(UnLock, self).promise_me(ContainerVector(*args))