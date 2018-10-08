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


def promise_me(module, *args):
    arg_vector = libMA.VectorPledge()
    for arg in args:
        arg_vector.append(arg)
    if isinstance(module, libMA.Module):
        return libMA.ModulePledge(module, arg_vector)
    elif isinstance(module, libMA.VolatileModule):
        return libMA.VolatileModulePledge(module, arg_vector)
    else:
        raise Exception("module must be an instance of Module or VolatileModule")

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
    # @brief Execute the implemented algorithm.
    # @details
    # Reimplemented from @ref Module::saveExecute.
    def execute(self, *args):
        self.__store_result = Module.execute(args)
        return self.__store_result

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
# @brief Contains sets of seeds ordered by SoC score
# @details
#
# @ingroup container
#
class SoCPriorityQueue(libMA.SoCPriorityQueue):
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
            for index in range(val.start or 0, val.stop or len(self), val.step
                               or 1):
                ret.append(super(NucSeq, self).__getitem__(index))
            return ret
        else:
            return super(NucSeq, self).__getitem__(val)


##
# @brief python wrapper for BinarySeeding
class BinarySeeding(libMA.BinarySeeding):
    def execute(self, *args):
        return super(BinarySeeding, self).execute(ContainerVector(list(args)))

    def promise_me(self, *args):
        return super(BinarySeeding, self).promise_me(
            ContainerVector(list(args)))


##
# @brief python wrapper for FileWriter
class FileWriter(libMA.FileWriter):
    def execute(self, *args):
        return super(FileWriter, self).execute(ContainerVector(list(args)))

    def promise_me(self, *args):
        return super(FileWriter, self).promise_me(ContainerVector(list(args)))


##
# @brief python wrapper for PairedReads
# class PairedReads(libMA.PairedReads):
#     def execute(self, *args):
#         return super(PairedReads, self).execute(ContainerVector(list(args)))
# 
#     def promise_me(self, *args):
#         return super(PairedReads, self).promise_me(ContainerVector(list(args)))


##
# @brief python wrapper for StripOfConsideration
class StripOfConsideration(libMA.StripOfConsideration):
    def execute(self, *args):
        return super(StripOfConsideration, self).execute(
            ContainerVector(list(args)))

    def promise_me(self, *args):
        return super(StripOfConsideration, self).promise_me(
            ContainerVector(list(args)))


##
# @brief python wrapper for Harmonization
class Harmonization(libMA.Harmonization):
    def execute(self, *args):
        return super(Harmonization, self).execute(
            ContainerVector(list(args)))

    def promise_me(self, *args):
        return super(Harmonization, self).promise_me(
            ContainerVector(list(args)))


##
# @brief python wrapper for FileReader
class FileReader(libMA.FileReader):
    def execute(self):
        return super(FileReader, self).execute(ContainerVector(Nil()))

    def promise_me(self, *args):
        return super(FileReader, self).promise_me(ContainerVector(list(args)))


##
# @brief python wrapper for PairedFileReader
# class PairedFileReader(libMA.PairedFileReader):
#     def execute(self):
#         return super(PairedFileReader, self).execute(ContainerVector(Nil()))
# 
#     def promise_me(self, *args):
#         return super(PairedFileReader, self).promise_me(
#             ContainerVector(list(args)))


##
# @brief python wrapper for LinearLineSweep
class NeedlemanWunsch(libMA.NeedlemanWunsch):
    def execute(self, *args):
        return super(NeedlemanWunsch, self).execute(
            ContainerVector(list(args)))

    def promise_me(self, *args):
        return super(NeedlemanWunsch, self).promise_me(
            ContainerVector(list(args)))


##
# @brief python wrapper for MappingQuality
class MappingQuality(libMA.MappingQuality):
    def execute(self, *args):
        return super(MappingQuality, self).execute(ContainerVector(list(args)))

    def promise_me(self, *args):
        return super(MappingQuality, self).promise_me(
            ContainerVector(list(args)))


##
# @brief The Lock Module.
# @ingroup module
#
class Lock(libMA.Lock):
    def execute(self, *args):
        return super(Lock, self).execute(ContainerVector(list(args)))

    def promise_me(self, *args):
        return super(Lock, self).promise_me(ContainerVector(list(args)))


##
# @brief The UnLock Module.
# @ingroup module
#
class UnLock(libMA.UnLock):
    def execute(self, *args):
        return super(UnLock, self).execute(ContainerVector(list(args)))

    def promise_me(self, *args):
        return super(UnLock, self).promise_me(ContainerVector(list(args)))


##
# @brief The OtherSeeding Module.
# @ingroup module
#
class OtherSeeding(libMA.OtherSeeding):
    def execute(self, *args):
        return super(OtherSeeding, self).execute(ContainerVector(list(args)))

    def promise_me(self, *args):
        return super(OtherSeeding, self).promise_me(
            ContainerVector(list(args)))


#def ksw(a_list, b_list, m=2, mm=4, o=4, e=2, o2=24, e2=1, w=0):
#    ret = []
#    for a, b in zip(a_list, b_list):
#        cigar = libMA.run_ksw(a, b, m, mm, o, e, o2, e2, w)
#        cigar_list = []
#        for x in cigar.split(",")[:-1]:
#            num = int(x[:-1])
#            sym = x[-1]
#            cigar_list.append((num, sym))
#        ret.append(cigar_list)
#    return ret