#include "container/pack.h"
using namespace libMABS;

void exportPack()
{
    boost::python::class_<
            Pack, 
            boost::noncopyable,
            boost::python::bases<Container>,
            std::shared_ptr<Pack>
        >(
                "Pack",
                "unpacked_size_single_strand: the size of one strand\n"
                "\n"
                "Holds a packed sequence.\n"
                "Usefull for long sequences.\n"
            )
        .def(
                "unpacked_size", 
                &Pack::uiUnpackedSizeForwardPlusReverse,
                "arg1: self\n"
                "returns: the length of the sequence within the pack.\n"
            )
        .def(
                "append", 
                &Pack::vAppendSequence_boost,
                "arg1: self\n"
                "arg2: a NucSeq\n"
                "returns: nil\n"
                "\n"
                "Appends seq at the end of the pack.\n"
            )
#if 0 //DEPRECATED
        .def(
                "append_fasta_file", 
                &Pack::vAppendFastaFile,
                "arg1: self\n"
                "arg2: the filename\n"
                "returns: nil\n"
                "\n"
                "Appends seq at the end of the pack.\n"
            )
#endif
        .def(
                "store", 
                &Pack::vStoreCollection,
                "arg1: self\n"
                "arg2: the folder and filename on disk\n"
                "returns: nil\n"
                "\n"
                "Stores this pack at the given location.\n"
            )
        .def(
                "exists", 
                &Pack::packExistsOnFileSystem,
                "arg1: self\n"
                "arg2: the folder and filename on disk\n"
                "returns: a bool indicating if the file exists\n"
                "\n"
                "Checks weather a pack exists at the given location.\n"
            )
        .staticmethod("exists")
        .def(
                "load", 
                &Pack::vLoadCollection,
                "arg1: self\n"
                "arg2: the folder and filename on disk\n"
                "returns: nil\n"
                "\n"
                "Loads a pack from the given location on disk.\n"
            )
        .def(
                "extract_from_to", 
                &Pack::vExtract,
                "arg1: self\n"
                "arg2: begin of extraction\n"
                "arg3: end of extraction\n"
                "returns: the extracted sequence as NucSeq\n"
                "\n"
                "Extracts a sequence from the pack.\n"
                "Indices are inclusive.\n"
            )
        .def(
                "extract_complete", 
                &Pack::vColletionAsNucSeq,
                "arg1: self\n"
                "returns: the extracted sequence as NucSeq\n"
                "\n"
                "Extracts the entire pack as sequence.\n"
            )
        .def(
                "extract_forward_strand", 
                &Pack::vColletionWithoutReverseStrandAsNucSeq,
                "arg1: self\n"
                "returns: the extracted sequence as NucSeq\n"
                "\n"
                "Extracts the forward strand of the pack as sequence.\n"
            )
        .def(
                "extract_reverse_strand", 
                &Pack::vColletionOnlyReverseStrandAsNucSeq,
                "arg1: self\n"
                "returns: the extracted sequence as NucSeq\n"
                "\n"
                "Extracts the reverse strand of the pack as sequence.\n"
            )
        .def(
                "is_bridging", 
                &Pack::bridgingSubsection_boost
            )
        .def(
                "start_of_sequence", 
                &Pack::startOfSequenceWithName
            )
        .def(
                "start_of_sequence_id", 
                &Pack::startOfSequenceWithId
            )
        .def_readonly(
                "unpacked_size_single_strand", 
                &Pack::uiUnpackedSizeForwardStrand
            )
        ;

    //tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible<
            std::shared_ptr<Pack>,
            std::shared_ptr<Container>
        >();
}//function
