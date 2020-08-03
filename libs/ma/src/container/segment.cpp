/**
 * @file segment.cpp
 * @author Markus Schmidt
 */
#include "ma/container/segment.h"
#include "ms/util/pybind11.h"
using namespace libMA;

#ifdef WITH_PYTHON

void exportSegment( libMS::SubmoduleOrganizer& xOrganizer )
{
    // export the SegmentVector class
    py::class_<Segment, libMS::Container, std::shared_ptr<Segment>>( xOrganizer.container( ), "Segment" )
        .def( "start", &Segment::start_boost1 )
        .def( "end", &Segment::end_boost1 );

    // export the SegmentVector class
    py::bind_vector_ext<SegmentVector, libMS::Container, std::shared_ptr<SegmentVector>>( xOrganizer.container( ),
                                                                                          "SegmentVector", "docstr" )
        .def( "extract_seeds", &SegmentVector::extractSeeds )
        .def( "num_seeds", &SegmentVector::numSeeds )
#if MEASURE_DURATIONS == ( 1 )
        .def_readwrite( "fExtraction", &SegmentVector::fExtraction )
        .def_readwrite( "fSorting", &SegmentVector::fSorting )
        .def_readwrite( "fLinesweep", &SegmentVector::fLinesweep )
#endif
        .def( "num_seeds_larger", &SegmentVector::numSeedsLarger );

    // tell boost python that pointers of these classes can be converted implicitly
    py::implicitly_convertible<SegmentVector, libMS::Container>( );

} // function
#endif