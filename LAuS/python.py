##
# @package LAuS.python
# @brief This package is merely for documentation purpose.
# @file python.py
# @brief This file is merely for documentation purpose.
# @details
# This file does not get loaded nor installed however doxygen does not know that.
# Thus we can provide documentation for classes that exist 
# (since they are generated in the C++ code),
# but are not present in for of a individual file.
# @author Markus Schmidt
#

raise ImportError("You imported a Module that is used purely for documentation.")

from libLAuS import *
from collections.abc import MutableSequence

##
# @brief contains the final output of the aligner.
# @details
# Holds a sparse vector like representation of one alignment.
# @note To create an instance of this class omit the trainling underscore.
# @ingroup container
#
class Alignment_(Alignment):

    ##
    # @brief Get the MatchType at the position index.
    # @details
    # Reimplemented from @ref Alignment::at.
    #
    def at(self, index):
        pass

    ##
    # @brief Append amount elements of the given math_type to the end of the alignment.
    # @details
    # Reimplemented from @ref Alignment::append_boost1.
    #
    def append(self, match_type, amount):
        pass

    ##
    # @brief Append one element of the given math_type to the end of the alignment.
    # @details
    # Reimplemented from @ref Alignment::append_boost2.
    #
    def append(self, match_type):
        pass

    ##
    # @brief The starting position of the alignment on the reference.
    # @details
    # Reimplemented from @ref Alignment::beginOnRef.
    #
    def begin_on_ref(self):
        pass

    ##
    # @brief The ending position of the alignment on the reference.
    # @details
    # Reimplemented from @ref Alignment::endOnRef.
    #
    def end_on_ref(self):
        pass

    ##
    # @brief The length of the alignment.
    # @details
    # Reimplemented from @ref Alignment::length.
    #
    def __len__(self):
        pass

    ##
    # @brief The length of the alignment.
    # @details
    # Reimplemented from @ref Alignment::length.
    #
    def length(self):
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
# @note To create an instance of this class omit the trainling underscore.
#
class MatchType_(MatchType):
    def __init__(self):
        pass

##
# @brief representation of a interval
#
class Interval_(Interval):
    def __init__(self, start, size):
        self.start = start
        self.size = size

##
# @brief Represents the pledge to deliver some container.
# @details
# Content may be provided by a @ref CppModule::promiseMe "module" or by calling @ref set.
# @note To create an instance of this class omit the trainling underscore.
# @ingroup container
#
class Pledge_(Pledge):
    ##
    # @brief Create a pledge with a module thats responsible to fullfill it.
    # @details
    # This is used the the @ref Module::promise_me "promise_me" function of Module. <br>
    # Reimplemented from @ref Pledge.
    #
    def __init__(self, module, container_type, in_pledge):
        pass

    ##
    # @brief Manually create a pledge.
    # @details
    # You are required to fullfill you pledge by calling @ref Pledge_::set "set". <br>
    # Reimplemented from @ref Pledge.
    #
    def __init__(self, container_type):
        pass

    ##
    # @brief Manually fullfill a pledge.
    # @details
    # Reimplemented from @ref Pledge::set.
    #
    def set(self, container):
        pass

    ##
    # @brief Trigger the computational graph to fullfill this pledge
    # @details
    # This call will trigger all containers to invalidate their content.
    # Then every part of the computational graph that is needed to fullfill 
    # this pledge will be computed. <br>
    # Reimplemented from @ref Pledge::next.
    #
    def next(self, container):
        pass

    ##
    # @brief Trigger the computational graph to fullfill this pledge
    # @details
    # Every part of the computational graph that is needed to fullfill
    # this pledge will be computed. <br>
    # Reimplemented from @ref Pledge::get.
    #
    def get(self, container):
        pass

##
# @brief Contains a suffix array.
# @details
# @note To create an instance of this class omit the trainling underscore.
# @ingroup container
#
class FMIndex_(FMIndex):
    ##
    # @brief Create a empty suffix array.
    # @details
    # Reimplemented from @ref FMIndex.
    #
    def __init__(self):
        pass

    ##
    # @brief Create a new suffix array using the given sequence.
    # @details
    # Reimplemented from @ref FMIndex.
    #
    def __init__(self, NucSeq):
        pass

    ##
    # @brief Create a new suffix array using the given sequence.
    # @details
    # Reimplemented from @ref FMIndex.
    #
    def __init__(self, Pack):
        pass

    ##
    # @brief Load a suffix array from disk.
    # @details
    # Reimplemented from @ref FMIndex::vLoadFMIndex.
    #
    def load(self, file_name):
        pass

    ##
    # @brief Check if there is a suffixarray named file_name on disk.
    # @details
    # Reimplemented from @ref FMIndex::packExistsOnFileSystem.
    #
    @staticmethod
    def exists(file_name):
        pass

    ##
    # @brief Store a suffix array on disk.
    # @details
    # Reimplemented from @ref FMIndex::vStoreFMIndex.
    #
    def store(self, file_name):
        pass

    ##
    # @brief Delivers the Position in the reference sequence that belongs to 
    # the position k in the BWT.
    # @details
    # Reimplemented from @ref FMIndex::bwt_sa.
    #
    def bwt_sa(self, k):
        pass

    ##
    # @brief Delivers the Position in the reference sequence that belongs to 
    # the position k in the BWT.
    # @details
    # Reimplemented from @ref FMIndex::bwt_2occ4.
    #
    def bwt_2occ4(self, k):
        pass

##
# @brief The Container for a SegmentVector
# @details
# A doubly linked list holding Segments.
# @note To create an instance of this class omit the trainling underscore.
# @ingroup container
#
class SegmentVector_(SegmentVector, MutableSequence):
    pass

##
# @brief The Container for a Seeds_
# @details
# is iterable
# @note To create an instance of this class omit the trainling underscore.
# @ingroup container
#
class Seeds_(Seeds, MutableSequence):
    pass

##
# @brief A interval on the query that contains one or multiple @ref Seed "seeds".
# @details
# @note To create an instance of this class omit the trainling underscore.
# @ingroup container
#
class Segment_(Segment, Interval_):
    ##
    # @brief Create a new Segment
    # @details
    # Reimplemented from @ref Segment.
    #
    def __init__(self, start, size, sa_interval):
        self.sa_interval = sa_interval

##
# @brief A packed version of a @ref NucSeq_ "NucSeq".
# @details
# @note To create an instance of this class omit the trainling underscore.
# @ingroup container
#
class Pack_(Pack):
    ##
    # @brief Create an empty pack.
    # @details
    # Reimplemented from @ref Pack.
    #
    def __init__(self):
        self.unpacked_size_single_strand = 0
        pass

    ##
    # @brief The length of the sequence within the pack
    # @details
    # Reimplemented from @ref Pack::uiUnpackedSizeForwardPlusReverse.
    #
    def unpacked_size(self):
        pass

    ##
    # @brief Append a @ref NucSeq_ "NucSeq".
    # @details
    # Reimplemented from @ref Pack::vAppendSequence_boost.
    #
    def append(self, nuc_seq):
        pass

    ##
    # @brief Stores this pack at the given location.
    # @details
    # Reimplemented from @ref Pack::vStoreCollection.
    #
    def store(self, file_name):
        pass

    ##
    # @brief Checks weather a pack exists at the given location.
    # @details
    # Reimplemented from @ref Pack::packExistsOnFileSystem.
    #
    @staticmethod
    def exists(file_name):
        pass

    ##
    # @brief Loads a pack from the given location on disk.
    # @details
    # Reimplemented from @ref Pack::vLoadCollection.
    #
    def load(self, file_name):
        pass

    ##
    # @brief Extracts a sequence from the pack.
    # @details
    # Indices are inclusive <br>
    # Reimplemented from @ref Pack::vExtract.
    #
    def extract_from_to(self, start, end):
        pass

    ##
    # @brief Extracts the entire pack as sequence.
    # @details
    # Reimplemented from @ref Pack::vColletionAsNucSeq.
    #
    def extract_complete(self):
        pass

    ##
    # @brief Extracts the forward strand entire pack as sequence.
    # @details
    # Reimplemented from @ref Pack::vColletionWithoutReverseStrandAsNucSeq.
    #
    def extract_forward_strand(self):
        pass

    ##
    # @brief Extracts the reverse strand entire pack as sequence.
    # @details
    # Reimplemented from @ref Pack::vColletionOnlyReverseStrandAsNucSeq.
    #
    def extract_reverse_strand(self):
        pass

##
# @brief A single seed.
# @details
# @note To create an instance of this class omit the trainling underscore.
# @ingroup container
#
class Seed_(Seed, Interval_):
    ##
    # @brief Create a new Seed.
    # @details
    # Reimplemented from @ref Seed.
    #
    def __init__(self, query_start, size, ref_start):
        self.ref_start = ref_start

##
# @brief A SAInterval.
# @details
# @note To create an instance of this class omit the trainling underscore.
# @ingroup container
#
class SAInterval_(SAInterval, Interval_):
    ##
    # @brief Create a new SAInterval.
    # @details
    # Reimplemented from @ref SAInterval.
    #
    def __init__(self, query_start, size):
        pass

##
# @brief A nucleotide sequence.
# @details
# @note To create an instance of this class omit the trainling underscore.
# @ingroup container
#
class NucSeq_(NucSeq):
    ##
    # @brief Create a nucleotide sequence out of the given string.
    # @details
    # Reimplemented from @ref NucSeq.
    #
    def __init__(self, string):
        pass

    ##
    # @brief Get the nucleotide at the given index.
    # @details
    # Reimplemented from @ref NucSeq::charAt.
    #
    def at(self, index):
        pass

    ##
    # @brief Get the nucleotide at the given index.
    # @details
    # Reimplemented from @ref NucSeq::charAt.
    #
    def __getitem__(self, index):
        pass

    ##
    # @brief Append the given sequence.
    # @details
    # Reimplemented from @ref NucSeq::vAppend_boost.
    #
    def append(self, string):
        pass

    ##
    # @brief Get the length of the sequence.
    # @details
    # Reimplemented from @ref NucSeq::length.
    #
    def length(self):
        pass

    ##
    # @brief Get the length of the sequence.
    # @details
    # Reimplemented from @ref NucSeq::length.
    #
    def __len__(self):
        pass

    ##
    # @brief Get sequence as a string.
    # @details
    # Reimplemented from @ref NucSeq::toString.
    #
    def __str__(self):
        pass

    ##
    # @brief Reverse the sequence.
    # @details
    # Reimplemented from @ref NucSeq::vReverse.
    #
    def reverse(self):
        pass

##
# @brief The StripOfConsideration Module.
# @details
# Used to quickly find areas with high density of @ref Seed "seeds".
# @note To create an instance of this class omit the trainling underscore.
# @ingroup module
#
class StripOfConsideration_(StripOfConsideration):
    ##
    # @brief Create a new Module.
    # @details
    # Reimplemented from @ref StripOfConsideration.
    #
    def __init__(self):
        ##
        # @brief The strip of consideration size.
        self.strip_size = 10000
        ##
        # @brief Maximum ambiguity for a seed to be considered.
        self.max_hits = 500
        ##
        # @brief skip seeds with too much ambiguity
        # @details
        # True: skip all seeds with to much ambiguity
        # False: use max_hits instances of the seeds with more ambiguity
        self.skip_long = True

    ##
    # @brief returns the @ref ContainerType "container types" segmentList, segmentList, nucSeq, packedNucSeq, fm_index.
    # @details
    # Reimplemented from @ref StripOfConsideration::getInputType.
    def get_input_type(self):
        pass

    ##
    # @brief returns the @ref ContainerType "container type" seedsVector.
    # @details
    # Reimplemented from @ref StripOfConsideration::getOutputType.
    def get_output_type(self):
        pass

#DEPRECATED
"""
##
# @brief The LineSweep Module.
# @details
# @note To create an instance of this class omit the trainling underscore.
# @ingroup module
#
class LineSweep_(LineSweep):
    ##
    # @brief Create a new Module.
    # @details
    # Reimplemented from @ref LineSweep.
    #
    def __init__(self):
        pass

    ##
    # @brief returns the @ref ContainerType "container types" nucSeq, packedNucSeq, seeds.
    # @details
    # Reimplemented from @ref LineSweep::getInputType.
    def get_input_type(self):
        pass

    ##
    # @brief returns the @ref ContainerType "container type" seeds.
    # @details
    # Reimplemented from @ref LineSweep::getOutputType.
    def get_output_type(self):
        pass
"""

##
# @brief The LinearLineSweep Module.
# @details
# @note To create an instance of this class omit the trainling underscore.
# @ingroup module
#
class LinearLineSweep_(LinearLineSweep):
    ##
    # @brief Create a new Module.
    # @details
    # Reimplemented from @ref LinearLineSweep.
    #
    def __init__(self):
        pass

    ##
    # @brief returns the @ref ContainerType "container types" nucSeq, packedNucSeq, seeds.
    # @details
    # Reimplemented from @ref LinearLineSweep::getInputType.
    # Expects:
    # - Seeds
    def get_input_type(self):
        pass

    ##
    # @brief returns the @ref ContainerType "container type" seeds.
    # @details
    # Reimplemented from @ref LinearLineSweep::getOutputType.
    # Type is: Seeds
    def get_output_type(self):
        pass

##
# @brief The Chaining Module.
# @details
# @note To create an instance of this class omit the trainling underscore.
# @ingroup module
#
class Chaining_(Chaining):
    ##
    # @brief Create a new Module.
    # @details
    # Reimplemented from @ref Chaining.
    #
    def __init__(self):
        pass

    ##
    # @brief returns the @ref ContainerType "container types" nucSeq, packedNucSeq, seeds.
    # @details
    # Reimplemented from @ref Chaining::getInputType.
    def get_input_type(self):
        pass

    ##
    # @brief returns the @ref ContainerType "container type" seeds.
    # @details
    # Reimplemented from @ref Chaining::getOutputType.
    def get_output_type(self):
        pass

##
# @brief The NeedlemanWunsch Module.
# @details
# @note To create an instance of this class omit the trainling underscore.
# @ingroup module
#
class NeedlemanWunsch_(NeedlemanWunsch):
    ##
    # @brief Create a new Module.
    # @details
    # Reimplemented from @ref NeedlemanWunsch.
    #
    def __init__(self):
        pass

    ##
    # @brief returns the @ref ContainerType "container types" seeds, nucSeq, packedNucSeq.
    # @details
    # Reimplemented from @ref NeedlemanWunsch::getInputType.
    def get_input_type(self):
        pass

    ##
    # @brief returns the @ref ContainerType "container type" alignment.
    # @details
    # Reimplemented from @ref NeedlemanWunsch::getOutputType.
    def get_output_type(self):
        pass

##
# @brief The GetAnchors Module.
# @details
# @note To create an instance of this class omit the trainling underscore.
# @ingroup module
#
class GetAnchors_(GetAnchors):
    ##
    # @brief Create a new Module.
    # @details
    # Reimplemented from @ref GetAnchors.
    #
    def __init__(self):
        pass

    ##
    # @brief returns the @ref ContainerType "container types" segmentList.
    # @details
    # Reimplemented from @ref GetAnchors::getInputType.
    def get_input_type(self):
        pass

    ##
    # @brief returns the @ref ContainerType "container type" segmentList.
    # @details
    # Reimplemented from @ref GetAnchors::getOutputType.
    def get_output_type(self):
        pass

##
# @brief The BinarySeeding Module.
# @details
# @note To create an instance of this class omit the trainling underscore.
# @ingroup module
#
class BinarySeeding_(BinarySeeding):
    ##
    # @brief Create a new Module.
    # @details
    # Reimplemented from @ref BinarySeeding.
    # if use_lr_extension is true then TODO:
    #
    def __init__(self, use_lr_extension):
        pass

    ##
    # @brief returns the @ref ContainerType "container types" fm_index, nucSeq.
    # @details
    # Reimplemented from @ref BinarySeeding::getInputType.
    #
    def get_input_type(self):
        pass

    ##
    # @brief returns the @ref ContainerType "container type" segmentList.
    # @details
    # Reimplemented from @ref BinarySeeding::getOutputType.
    #
    def get_output_type(self):
        pass

##
# @brief The ContainerVector Module.
# @details
# Holds multiple containers.
# @note To create an instance of this class omit the trainling underscore.
# @ingroup container
#
class ContainerVector_(ContainerVector, MutableSequence):
    pass

