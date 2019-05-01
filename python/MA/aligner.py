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

from libMA import Alignment
from libMA import MatchType
from libMA import Pledge
from libMA import SoCPriorityQueue
from libMA import SegmentVector
from libMA import Seeds
from libMA import Segment
from libMA import Pack
from libMA import Seed
from libMA import SAInterval
from libMA import FMIndex
from libMA import Module
from libMA import VectorPledge
from libMA import VolatileModule
from libMA import configureAccurate
from libMA import configureFast
from libMA import NucSeq
from libMA import SV_DB
from libMA import SvJumpInserter
from libMA import SvCallInserter
from libMA import SvJump
from libMA import SvCall
from libMA import VectorPledge
from libMA import ReadInserter
from libMA import ParameterSetManager
from libMA import ExecutionContext
from libMA import ContainerVectorNucSeq
from libMA import SortedSvJumpFromSql
from libMA import AllNucSeqFromSql
from libMA import SvCallerRunsFromDb
from libMA import SvCallsFromDb
import libMA


def promise_me(module, *args):
    arg_vector = libMA.VectorPledge()
    for arg in args:
        arg_vector.append(arg)
    if isinstance(module, Module):
        return libMA.ModulePledge(module, arg_vector)
    elif isinstance(module, VolatileModule):
        return libMA.VolatileModulePledge(module, arg_vector)
    else:
        raise Exception(
            "module must be an instance of Module or VolatileModule")


##
# @brief A nucleotide sequence.
# @details
#
# @ingroup container
#
# class NucSeq(libMA.NucSeq):
#     ##
#     # @brief enable slicing a NucSeq
#     def __getitem__(self, val):
#         if isinstance(val, slice):
#             ret = NucSeq()
#             for index in range(val.start or 0, val.stop or len(self), val.step
#                                or 1):
#                 ret.append(super(NucSeq, self).__getitem__(index))
#             return ret
#         else:
#             return super(NucSeq, self).__getitem__(val)

##
# @brief python wrapper for BinarySeeding
class BinarySeeding(libMA.BinarySeeding):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(BinarySeeding, self).execute(vec)


##
# @brief python wrapper for FileWriter
class FileWriter(libMA.FileWriter):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(FileWriter, self).execute(vec)


##
# @brief python wrapper for PairedFileWriter
class PairedFileWriter(libMA.PairedFileWriter):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(PairedFileWriter, self).execute(vec)


if hasattr(libMA, "DBWriter"):
    ##
    # @brief python wrapper for DbWriter
    # @details
    # This requires a switch in order to be compiled
    class DbWriter(libMA.DbWriter): #pylint: disable=I1101
        def execute(self, *args):
            vec = libMA.ContainerVector()
            for arg in args:
                vec.append(arg)
            return super(DbWriter, self).execute(vec)
else:
    ##
    # @brief python wrapper for DbWriter
    class DbWriter:
        def __init__(self):
            raise NotImplementedError(
                "libMA was complied without the WITH_POSTGRES" +
                " switch, so this class is unavailable")


if hasattr(libMA, "PairedDbWriter"):
    ##
    # @brief python wrapper for PairedDbWriter
    # This requires a switch in order to be compiled
    class PairedDbWriter(libMA.PairedDbWriter): #pylint: disable=I1101
        def execute(self, *args):
            vec = libMA.ContainerVector()
            for arg in args:
                vec.append(arg)
            return super(PairedDbWriter, self).execute(vec)
else:
    ##
    # @brief python wrapper for PairedDbWriter
    class PairedDbWriter:
        def __init__(self):
            raise NotImplementedError(
                "libMA was complied without the WITH_POSTGRES" +
                " switch, so this class is unavailable")


##
# @brief python wrapper for PairedReads
class PairedReads(libMA.PairedReads):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(PairedReads, self).execute(vec)


##
# @brief python wrapper for StripOfConsideration
class StripOfConsideration(libMA.StripOfConsideration):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(StripOfConsideration, self).execute(vec)


##
# @brief python wrapper for Harmonization
class Harmonization(libMA.Harmonization):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(Harmonization, self).execute(vec)


##
# @brief python wrapper for HarmonizationSingle
class HarmonizationSingle(libMA.HarmonizationSingle):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(HarmonizationSingle, self).execute(vec)


##
# @brief python wrapper for SeedLumping
class SeedLumping(libMA.SeedLumping):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(SeedLumping, self).execute(vec)


##
# @brief python wrapper for FileReader
class FileReader(libMA.FileReader):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(FileReader, self).execute(vec)


##
# @brief python wrapper for PairedFileReader
class PairedFileReader(libMA.PairedFileReader):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(PairedFileReader, self).execute(vec)


##
# @brief python wrapper for LinearLineSweep
class NeedlemanWunsch(libMA.NeedlemanWunsch):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(NeedlemanWunsch, self).execute(vec)


##
# @brief python wrapper for MappingQuality
class MappingQuality(libMA.MappingQuality):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(MappingQuality, self).execute(vec)


##
# @brief The Lock Module.
# @ingroup module
#
class Lock(libMA.Lock):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(Lock, self).execute(vec)


##
# @brief The UnLock Module.
# @ingroup module
#
class UnLock(libMA.UnLock):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(UnLock, self).execute(vec)


##
# @brief The UnLock Module.
# @ingroup module
#
class GetFirstQuery(libMA.GetFirstQuery):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(GetFirstQuery, self).execute(vec)


##
# @brief The UnLock Module.
# @ingroup module
#
class GetSecondQuery(libMA.GetSecondQuery):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(GetSecondQuery, self).execute(vec)


##
# @brief The OtherSeeding Module.
# @ingroup module
#
class OtherSeeding(libMA.OtherSeeding):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(OtherSeeding, self).execute(vec)


# @brief The SmallInversions Module.
# @ingroup module
#
class SmallInversions(libMA.SmallInversions):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(SmallInversions, self).execute(vec)


##
# @brief The NucSeqFromSql Module.
# @ingroup module
#
class NucSeqFromSql(libMA.NucSeqFromSql):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(NucSeqFromSql, self).execute(vec)


##
# @brief The AlignmentCollector Module.
# @ingroup module
#
class AlignmentCollector(libMA.AlignmentCollector):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(NucSeqFromSql, self).execute(vec)


##
# @brief The PairedNucSeqFromSql Module.
# @ingroup module
#
class PairedNucSeqFromSql(libMA.PairedNucSeqFromSql):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(PairedNucSeqFromSql, self).execute(vec)


##
# @brief The NucSeqSplitter Module.
# @ingroup module
#
class NucSeqSplitter(libMA.NucSeqSplitter):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(NucSeqSplitter, self).execute(vec)


##
# @brief The FillSeedSet Module.
# @ingroup module
#
class FillSeedSet(libMA.FillSeedSet):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(FillSeedSet, self).execute(vec)


##
# @brief The SvDbInserter Module.
# @ingroup module
#
class SvDbInserter(libMA.SvDbInserter):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(SvDbInserter, self).execute(vec)


##
# @brief The SvJumpsFromSeeds Module.
# @ingroup module
#
class SvJumpsFromSeeds(libMA.SvJumpsFromSeeds):
    def execute(self, *args):
        vec = libMA.ContainerVector()
        for arg in args:
            vec.append(arg)
        return super(SvJumpsFromSeeds, self).execute(vec)


##
# @brief convert bytes to a NucSeq
# @details
# Usefull for converting reads stored as blob data in sqlite3 to NucSeq objects.
#
def nuc_seq_from_bytes(blob):
    converter = libMA.NucSeqSql()
    converter.fromBlob(blob)
    return converter.seq
