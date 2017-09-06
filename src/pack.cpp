#include "pack.h"

void exportPack()
{
    boost::python::class_<
            BWACompatiblePackedNucleotideSequencesCollection, 
            boost::noncopyable,
            boost::python::bases<Container>,
            std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection>
        >(
                "BWAPack",
                "class: BWAPack\n"
                "   unpacked_size_single_strand: the size of one strand\n"
                "\n"
                "   Holds a packed sequence.\n"
                "   Usefull for long sequences.\n"
            )
        .def(
                "unpacked_size", 
                &BWACompatiblePackedNucleotideSequencesCollection::uiUnpackedSizeForwardPlusReverse,
                "method: unpacked_size\n"
                "   returns: the length of the sequence within the pack.\n"
            )
        .def(
                "append", 
                &BWACompatiblePackedNucleotideSequencesCollection::vAppendSequence_boost,
                "method: append(seq)\n"
                "   seq: a NucleotideSequence\n"
                "   returns: nil\n"
                "\n"
                "Appends seq at the end of the pack.\n"
            )
        .def(
                "store", 
                &BWACompatiblePackedNucleotideSequencesCollection::vStoreCollection,
                "method: store(file_prefix)\n"
                "   file_prefix: the folder and filename on disk\n"
                "   retuns: nil\n"
                "\n"
                "Stores this pack at the given loacation.\n"
            )
        .def(
                "exists", 
                &BWACompatiblePackedNucleotideSequencesCollection::packExistsOnFileSystem,
                "method: exists(file_prefix) static\n"
                "   file_prefix: the folder and filename on disk\n"
                "   retuns: a bool indicating if the file exists\n"
                "\n"
                "Checks weather a pack exists at the given location.\n"
            )
        .staticmethod("exists")
        .def(
                "load", 
                &BWACompatiblePackedNucleotideSequencesCollection::vLoadCollection,
                "method: load(file_prefix)\n"
                "   file_prefix: the folder and filename on disk\n"
                "   retuns: nil\n"
                "\n"
                "Loads a pack from the given location on disk.\n"
            )
        .def(
                "extract_from_to", 
                &BWACompatiblePackedNucleotideSequencesCollection::vExtract,
                "method: extract_from_to(from, to)\n"
                "   from: begin of extraction\n"
                "   to: end of extraction\n"
                "   retuns: the extracted sequence as NucleotideSequence\n"
                "\n"
                "Extracts a sequence from the pack.\n"
                "Indices are inclusive.\n"
            )
        .def(
                "extract_complete", 
                &BWACompatiblePackedNucleotideSequencesCollection::vColletionAsNucleotideSequence,
                "method: extract_complete()\n"
                "   retuns: the extracted sequence as NucleotideSequence\n"
                "\n"
                "Extracts the entire pack as sequence.\n"
            )
        .def(
                "extract_forward_strand", 
                &BWACompatiblePackedNucleotideSequencesCollection::vColletionWithoutReverseStrandAsNucleotideSequence,
                "method: extract_forward_strand()\n"
                "   retuns: the extracted sequence as NucleotideSequence\n"
                "\n"
                "Extracts the forward strand of the pack as sequence.\n"
            )
        .def(
                "extract_reverse_strand", 
                &BWACompatiblePackedNucleotideSequencesCollection::vColletionOnlyReverseStrandAsNucleotideSequence,
                "method: extract_reverse_strand()\n"
                "   retuns: the extracted sequence as NucleotideSequence\n"
                "\n"
                "Extracts the reverse strand of the pack as sequence.\n"
            )
        .def_readonly(
                "unpacked_size_single_strand", 
                &BWACompatiblePackedNucleotideSequencesCollection::uiUnpackedSizeForwardStrand
            )
        ;

	//tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible<
            std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection>,
            std::shared_ptr<Container>
        >();
}//function
