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
        .def("append", &BWACompatiblePackedNucleotideSequencesCollection::vAppendSequenceWrapper)
        .def("store", &BWACompatiblePackedNucleotideSequencesCollection::vStoreCollectionWrapper)
        .def("exists", &BWACompatiblePackedNucleotideSequencesCollection::packExistsOnFileSystemWrapper)
        .staticmethod("exits")
        .def("load", &BWACompatiblePackedNucleotideSequencesCollection::vLoadCollectionWrapper)
        .def("extractFromTo", &BWACompatiblePackedNucleotideSequencesCollection::vExtractSubsectionWrapper)
        .def("extractComplete", &BWACompatiblePackedNucleotideSequencesCollection::vColletionAsNucleotideSequenceWrapper)
        .def("extractForwardStrand", &BWACompatiblePackedNucleotideSequencesCollection::vColletionWithoutReverseStrandAsNucleotideSequenceWrapper)
        .def("extractReverseStrand", &BWACompatiblePackedNucleotideSequencesCollection::vColletionOnlyReverseStrandAsNucleotideSequenceWrapper)
        .def("unpackedSizeSingleStrand", &BWACompatiblePackedNucleotideSequencesCollection::uiStartOfReverseStrand);

        
	//tell boost python that it's possible to convert shared pointers with these classes
    boost::python::implicitly_convertible<std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection>,std::shared_ptr<Container>>();
}//function
