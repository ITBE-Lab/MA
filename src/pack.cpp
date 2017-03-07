#include "pack.h"

void exportPack()
{
    boost::python::class_<
            PackContainer, 
            boost::python::bases<Container>,
            std::shared_ptr<PackContainer>
        >("BWAPack")
        .def("unpackedSize", &PackContainer::getUnpackedSize)
        .def("append", &PackContainer::vAppendSequence)
        .def("store", &PackContainer::vStoreCollection)
        .def("exists", &PackContainer::packExistsOnFileSystem)
        .staticmethod("exists")
        .def("load", &PackContainer::vLoadCollection)
        .def("extractFromTo", &PackContainer::vExtractSubsection)
        .def("extractComplete", &PackContainer::vColletionAsNucleotideSequence)
        .def("extractForwardStrand", &PackContainer::vColletionWithoutReverseStrandAsNucleotideSequence)
        .def("extractReverseStrand", &PackContainer::vColletionOnlyReverseStrandAsNucleotideSequence)
        .def("unpackedSizeSingleStrand", &PackContainer::uiStartOfReverseStrand);

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< std::shared_ptr<PackContainer>, std::shared_ptr<Container> >(); 
}//function
