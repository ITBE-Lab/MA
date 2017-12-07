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

from libLAuS import ""

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
class FMIndex_(FM_Index):
    ##
    # @brief Create a empty suffix array.
    # @details
    # Reimplemented from @ref FM_Index.
    #
    def __init__(self):
        pass

    ##
    # @brief Create a new suffix array using the given sequence.
    # @details
    # Reimplemented from @ref FM_Index.
    #
    def __init__(self, NucSeq):
        pass

    ##
    # @brief Create a new suffix array using the given sequence.
    # @details
    # Reimplemented from @ref FM_Index.
    #
    def __init__(self, Pack):
        pass

    ##
    # @brief Load a suffix array from disk.
    # @details
    # Reimplemented from @ref FM_Index::vLoadFM_Index.
    #
    def load(self, file_name):
        pass

    ##
    # @brief Check if there is a suffixarray named file_name on disk.
    # @details
    # Reimplemented from @ref FM_Index::packExistsOnFileSystem.
    #
    @staticmethod
    def exists(file_name):
        pass

    ##
    # @brief Store a suffix array on disk.
    # @details
    # Reimplemented from @ref FM_Index::vStoreFM_Index.
    #
    def store(self, file_name):
        pass

    ##
    # @brief Delivers the Position in the reference sequence that belongs to 
    # the position k in the BWT.
    # @details
    # Reimplemented from @ref FM_Index::bwt_sa.
    #
    def bwt_sa(self, k):
        pass

    ##
    # @brief Delivers the Position in the reference sequence that belongs to 
    # the position k in the BWT.
    # @details
    # Reimplemented from @ref FM_Index::bwt_2occ4.
    #
    def bwt_2occ4(self, k):
        pass

##
# @brief The Container for a SegmentList
# @details
# A doubly linked list holding Segments.
# @note To create an instance of this class omit the trainling underscore.
# @ingroup container
#
class SegmentList_(SegmentTree):
    ##
    # @brief Create a empty list.
    # @details
    # Reimplemented from @ref SegmentTree.
    #
    def __init__(self):
        pass

    ##
    # @brief Create a list containing one Segment of the size query_length.
    # @details
    # Reimplemented from @ref SegmentTree.
    #
    def __init__(self, query_length):
        pass

    ##
    # @brief An SegmentListIterator_ pointing to the first element of the list.
    # @details
    # Reimplemented from @ref SegmentTree::begin.
    #
    def __iter__(self):
        pass

    ##
    # @brief Insert a segment before the position of the given iterator.
    # @details
    # Invalidates the iterator.
    # Replace it with the returned one. <br>
    # Reimplemented from @ref SegmentTree::insertBefore_boost.
    # @returns a new SegmentListIterator_.
    #
    def insert_before(self, segment, iterator):
        pass

    ##
    # @brief Insert a segment after the position of the given iterator.
    # @details
    # Invalidates the iterator.
    # Replace it with the returned one. <br>
    # Reimplemented from @ref SegmentTree::insertAfter_boost.
    # @returns a new SegmentListIterator_.
    #
    def insert_after(self, segment, iterator):
        pass

    ##
    # @brief Insert a segment at the end of the list.
    # @details
    # Reimplemented from @ref SegmentTree::push_back.
    #
    def push_back(self, segment):
        pass

    ##
    # @brief Insert a segment at the front of the list.
    # @details
    # Reimplemented from @ref SegmentTree::push_front.
    #
    def push_front(self, segment):
        pass

    ##
    # @brief Removes the element the SegmentListIterator_ points to.
    # @details
    # Invalidates the iterator. <br>
    # Reimplemented from @ref SegmentTree::removeNode_boost.
    #
    def remove_node(self, segment, iterator):
        pass

    ##
    # @brief Returns all Seeds.
    # @details
    # Reimplemented from @ref SegmentTree::getSeeds.
    #
    def get_seeds(self, fm_index):
        pass

    ##
    # @brief Returns the number of Seeds.
    # @details
    # Reimplemented from @ref SegmentTree::numSeeds.
    #
    def num_seeds(self):
        pass

##
# @brief Iterator for the doubly linked SegmentList.
# @details
# @note To create an instance of this class omit the trainling underscore.
#
class SegmentListIterator_(SegmentTree::Iterator):
    ##
    # @brief Returns the current element and increments the iterator.
    # @details
    # Reimplemented from @ref SegmentTree::Iterator::next_boost.
    #
    def __next__(self):
        pass

##
# @brief A interval on the query that contains one or multiple @ref Seed "seeds".
# @details
# @note To create an instance of this class omit the trainling underscore.
# @ingroup container
#
class Segment_(SegmentTreeInterval):
    ##
    # @brief Create a new segment.
    # @details
    # Reimplemented from @ref SegmentTreeInterval.
    #
    def __init__(self, start, size):
        pass

    ##
    # @brief Get the start of the segment.
    # @details
    # Reimplemented from @ref SegmentTreeInterval::start_boost1.
    #
    def start(self):
        pass

    ##
    # @brief Set the start of the segment.
    # @details
    # Reimplemented from @ref SegmentTreeInterval::start_boost2.
    #
    def start(self, new):
        pass

    ##
    # @brief Get the end of the segment.
    # @details
    # Reimplemented from @ref SegmentTreeInterval::end_boost1.
    #
    def end(self):
        pass

    ##
    # @brief Set the end of the segment.
    # @details
    # Reimplemented from @ref SegmentTreeInterval::end_boost2.
    #
    def end(self, new):
        pass

"""
    TODO:
    ##
    # @brief Record a bwt interval that matches this segment.
    # @details
    #
    # Reimplemented from @ref SegmentTreeInterval::push_back.
    #
    def push_back_bwt(start):
        pass
"""

    ##
    # @brief Extracts all seeds from the segment.
    # @details
    # Reimplemented from @ref SegmentTreeInterval::getSeeds.
    #
    def get_seeds(self, fm_index):
        pass

    ##
    # @brief Get the length of the segment.
    # @details
    # Reimplemented from @ref SegmentTreeInterval::size_boost1.
    #
    def size(self):
        pass

    ##
    # @brief Set the length of the segment.
    # @details
    # Reimplemented from @ref SegmentTreeInterval::size_boost2.
    #
    def size(self, new):
        pass

    ##
    # @brief Set the start and end of the segment.
    # @details
    # Reimplemented from @ref SegmentTreeInterval::set.
    #
    def set(self, start, end):
        pass

##
# @brief A packed version of a @ref NucSeq_ "NucSeq".
# @details
# @note To create an instance of this class omit the trainling underscore.
# @ingroup container
#
class Pack_(BWACompatiblePackedNucleotideSequencesCollection):
    ##
    # @brief Create an empty pack.
    # @details
    # Reimplemented from @ref BWACompatiblePackedNucleotideSequencesCollection.
    #
    def __init__(self):
        self.unpacked_size_single_strand = 0
        pass

    ##
    # @brief The length of the sequence within the pack
    # @details
    # Reimplemented from @ref BWACompatiblePackedNucleotideSequencesCollection::uiUnpackedSizeForwardPlusReverse.
    #
    def unpacked_size(self):
        pass

    ##
    # @brief Append a @ref NucSeq_ "NucSeq".
    # @details
    # Reimplemented from @ref BWACompatiblePackedNucleotideSequencesCollection::vAppendSequence_boost.
    #
    def append(self, nuc_seq):
        pass

    ##
    # @brief Stores this pack at the given location.
    # @details
    # Reimplemented from @ref BWACompatiblePackedNucleotideSequencesCollection::vStoreCollection.
    #
    def store(self, file_name):
        pass

    ##
    # @brief Checks weather a pack exists at the given location.
    # @details
    # Reimplemented from @ref BWACompatiblePackedNucleotideSequencesCollection::packExistsOnFileSystem.
    #
    @staticmethod
    def exists(file_name):
        pass

    ##
    # @brief Loads a pack from the given location on disk.
    # @details
    # Reimplemented from @ref BWACompatiblePackedNucleotideSequencesCollection::vLoadCollection.
    #
    def load(self, file_name):
        pass

    ##
    # @brief Extracts a sequence from the pack.
    # @details
    # Indices are inclusive <br>
    # Reimplemented from @ref BWACompatiblePackedNucleotideSequencesCollection::vExtract.
    #
    def extract_from_to(self, start, end):
        pass

    ##
    # @brief Extracts the entire pack as sequence.
    # @details
    # Reimplemented from @ref BWACompatiblePackedNucleotideSequencesCollection::vColletionAsNucleotideSequence.
    #
    def extract_complete(self):
        pass

    ##
    # @brief Extracts the forward strand entire pack as sequence.
    # @details
    # Reimplemented from @ref BWACompatiblePackedNucleotideSequencesCollection::vColletionWithoutReverseStrandAsNucleotideSequence.
    #
    def extract_forward_strand(self):
        pass

    ##
    # @brief Extracts the reverse strand entire pack as sequence.
    # @details
    # Reimplemented from @ref BWACompatiblePackedNucleotideSequencesCollection::vColletionOnlyReverseStrandAsNucleotideSequence.
    #
    def extract_reverse_strand(self):
        pass

##
# @brief A single seed.
# @details
# @note To create an instance of this class omit the trainling underscore.
# @ingroup container
#
class Seed_(Seed):
    ##
    # @brief The start position on the query.
    # @details
    # Reimplemented from @ref Seed::start_boost1.
    #
    def start(self):
        pass

    ##
    # @brief The end position on the query.
    # @details
    # Reimplemented from @ref Seed::end_boost1.
    #
    def end(self):
        pass

    ##
    # @brief The start position on the reference.
    # @details
    # Reimplemented from @ref Seed::start_ref.
    #
    def start_ref(self):
        pass

    ##
    # @brief The end position on the reference.
    # @details
    # Reimplemented from @ref Seed::end_ref.
    #
    def end_ref(self):
        pass

    ##
    # @brief The size of the seed.
    # @details
    # Reimplemented from @ref Seed::size_boost1.
    #
    def size(self):
        pass

##
# @brief A nucleotide sequence.
# @details
# @note To create an instance of this class omit the trainling underscore.
# @ingroup container
#
class NucSeq_(NucleotideSequence):
    ##
    # @brief Create a nucleotide sequence out of the given string.
    # @details
    # Reimplemented from @ref NucleotideSequence.
    #
    def __init__(self, string):
        pass

    ##
    # @brief Get the nucleotide at the given index.
    # @details
    # Reimplemented from @ref NucleotideSequence::charAt.
    #
    def at(self, index):
        pass

    ##
    # @brief Get the nucleotide at the given index.
    # @details
    # Reimplemented from @ref NucleotideSequence::charAt.
    #
    def __getitem__(self, index):
        pass

    ##
    # @brief Append the given sequence.
    # @details
    # Reimplemented from @ref NucleotideSequence::vAppend_boost.
    #
    def append(self, string):
        pass

    ##
    # @brief Get the length of the sequence.
    # @details
    # Reimplemented from @ref NucleotideSequence::length.
    #
    def length(self):
        pass

    ##
    # @brief Get the length of the sequence.
    # @details
    # Reimplemented from @ref NucleotideSequence::length.
    #
    def __len__(self):
        pass

    ##
    # @brief Get sequence as a string.
    # @details
    # Reimplemented from @ref NucleotideSequence::toString.
    #
    def __str__(self):
        pass

    ##
    # @brief Reverse the sequence.
    # @details
    # Reimplemented from @ref NucleotideSequence::vReverse.
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
class Chaining(Chaining_):
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
# @brief The NlongestIntervalsAsAnchors Module.
# @details
# @note To create an instance of this class omit the trainling underscore.
# @ingroup module
#
class NlongestIntervalsAsAnchors_(NlongestIntervalsAsAnchors):
    ##
    # @brief Create a new Module.
    # @details
    # Reimplemented from @ref NlongestIntervalsAsAnchors.
    #
    def __init__(self):
        pass

    ##
    # @brief returns the @ref ContainerType "container types" segmentList.
    # @details
    # Reimplemented from @ref NlongestIntervalsAsAnchors::getInputType.
    def get_input_type(self):
        pass

    ##
    # @brief returns the @ref ContainerType "container type" segmentList.
    # @details
    # Reimplemented from @ref NlongestIntervalsAsAnchors::getOutputType.
    def get_output_type(self):
        pass

##
# @brief The LongestNonEnclosedSegments Module.
# @details
# @note To create an instance of this class omit the trainling underscore.
# @ingroup module
#
class LongestNonEnclosedSegments_(LongestNonEnclosedSegments):
    ##
    # @brief Create a new Module.
    # @details
    # Reimplemented from @ref LongestNonEnclosedSegments.
    #
    def __init__(self):
        pass

    ##
    # @brief returns the @ref ContainerType "container types" fm_index, nucSeq.
    # @details
    # Reimplemented from @ref LongestNonEnclosedSegments::getInputType.
    #
    def get_input_type(self):
        pass

    ##
    # @brief returns the @ref ContainerType "container type" segmentList.
    # @details
    # Reimplemented from @ref LongestNonEnclosedSegments::getOutputType.
    #
    def get_output_type(self):
        pass

##
# @brief The LongestLRSegments Module.
# @details
# @note To create an instance of this class omit the trainling underscore.
# @ingroup module
#
class LongestLRSegments_(LongestLRSegments):
    ##
    # @brief Create a new Module.
    # @details
    # Reimplemented from @ref LongestLRSegments.
    #
    def __init__(self):
        pass

    ##
    # @brief returns the @ref ContainerType "container types" fm_index, nucSeq.
    # @details
    # Reimplemented from @ref LongestLRSegments::getInputType.
    #
    def get_input_type(self):
        pass

    ##
    # @brief returns the @ref ContainerType "container type" segmentList.
    # @details
    # Reimplemented from @ref LongestLRSegments::getOutputType.
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
class ContainerVector_(ContainerVector):
    ##
    # @brief Create a new Module.
    # @details
    # Reimplemented from @ref ContainerVector.
    #
    def __init__(self):
        pass

    ##
    # @brief Is iterable.
    #
    def __iter__(self):
        pass

    ##
    # @brief Is iterable.
    #
    def __len__(self):
        pass

    ##
    # @brief Is iterable.
    #
    def __at__(self, index):
        pass
