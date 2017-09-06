#include "pack.h"

void exportPack()
{
    boost::python::class_<
            BWACompatiblePackedNucleotideSequencesCollection, 
            boost::noncopyable,
            boost::python::bases<Container>,
            std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection>
        >("BWAPack")
        .def("unpackedSize", &BWACompatiblePackedNucleotideSequencesCollection::uiUnpackedSizeForwardPlusReverse)
        .def("append", &BWACompatiblePackedNucleotideSequencesCollection::vAppendSequence)
        .def("store", &BWACompatiblePackedNucleotideSequencesCollection::vStoreCollection)
        .def("exists", &BWACompatiblePackedNucleotideSequencesCollection::packExistsOnFileSystem)
        .staticmethod("exists")
        .def("load", &BWACompatiblePackedNucleotideSequencesCollection::vLoadCollection)
        .def("extractFromTo", &BWACompatiblePackedNucleotideSequencesCollection::vExtract)
        .def("extractComplete", &BWACompatiblePackedNucleotideSequencesCollection::vColletionAsNucleotideSequence)
        .def("extractForwardStrand", &BWACompatiblePackedNucleotideSequencesCollection::vColletionWithoutReverseStrandAsNucleotideSequence)
        .def("extractReverseStrand", &BWACompatiblePackedNucleotideSequencesCollection::vColletionOnlyReverseStrandAsNucleotideSequence)
        .def_readwrite("unpackedSizeSingleStrand", &BWACompatiblePackedNucleotideSequencesCollection::uiUnpackedSizeForwardStrand)
        ;

	//tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible<
            std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection>,
            std::shared_ptr<Container>
        >();
}//function
