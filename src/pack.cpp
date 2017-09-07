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
                "   unpacked_size_single_strand: the size of one strand\n"
                "\n"
                "Holds a packed sequence.\n"
                "Usefull for long sequences.\n"
            )
        .def(
                "unpacked_size", 
                &BWACompatiblePackedNucleotideSequencesCollection::uiUnpackedSizeForwardPlusReverse,
                "   arg1: self\n"
                "   returns: the length of the sequence within the pack.\n"
            )
        .def(
                "append", 
                &BWACompatiblePackedNucleotideSequencesCollection::vAppendSequence_boost,
                "   arg1: self\n"
                "   arg2: a NucleotideSequence\n"
                "   returns: nil\n"
                "\n"
                "Appends seq at the end of the pack.\n"
            )
        .def(
                "store", 
                &BWACompatiblePackedNucleotideSequencesCollection::vStoreCollection,
                "   arg1: self\n"
                "   arg2: the folder and filename on disk\n"
                "   returns: nil\n"
                "\n"
                "Stores this pack at the given loacation.\n"
            )
        .def(
                "exists", 
                &BWACompatiblePackedNucleotideSequencesCollection::packExistsOnFileSystem,
                "   arg1: self\n"
                "   arg2: the folder and filename on disk\n"
                "   returns: a bool indicating if the file exists\n"
                "\n"
                "Checks weather a pack exists at the given location.\n"
            )
        .staticmethod("exists")
        .def(
                "load", 
                &BWACompatiblePackedNucleotideSequencesCollection::vLoadCollection,
                "   arg1: self\n"
                "   arg2: the folder and filename on disk\n"
                "   returns: nil\n"
                "\n"
                "Loads a pack from the given location on disk.\n"
            )
        .def(
                "extract_from_to", 
                &BWACompatiblePackedNucleotideSequencesCollection::vExtract,
                "   arg1: self\n"
                "   arg2: begin of extraction\n"
                "   arg3: end of extraction\n"
                "   returns: the extracted sequence as NucleotideSequence\n"
                "\n"
                "Extracts a sequence from the pack.\n"
                "Indices are inclusive.\n"
            )
        .def(
                "extract_complete", 
                &BWACompatiblePackedNucleotideSequencesCollection::vColletionAsNucleotideSequence,
                "   arg1: self\n"
                "   returns: the extracted sequence as NucleotideSequence\n"
                "\n"
                "Extracts the entire pack as sequence.\n"
            )
        .def(
                "extract_forward_strand", 
                &BWACompatiblePackedNucleotideSequencesCollection::vColletionWithoutReverseStrandAsNucleotideSequence,
                "   arg1: self\n"
                "   returns: the extracted sequence as NucleotideSequence\n"
                "\n"
                "Extracts the forward strand of the pack as sequence.\n"
            )
        .def(
                "extract_reverse_strand", 
                &BWACompatiblePackedNucleotideSequencesCollection::vColletionOnlyReverseStrandAsNucleotideSequence,
                "   arg1: self\n"
                "   returns: the extracted sequence as NucleotideSequence\n"
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
