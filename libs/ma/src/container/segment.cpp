/**
 * @file segment.cpp
 * @author Markus Schmidt
 */
#include "container/segment.h"
#include "util/pybind11.h"
using namespace libMA;

#ifdef WITH_PYTHON

void exportSegment( py::module& rxPyModuleId )
{
    // export the SegmentVector class
    py::class_<Segment, Container, std::shared_ptr<Segment>>( rxPyModuleId, "Segment" )
        .def( "start", &Segment::start_boost1 )
        .def( "end", &Segment::end_boost1 );

    // export the SegmentVector class
    py::bind_vector_ext<SegmentVector, Container, std::shared_ptr<SegmentVector>>( rxPyModuleId, "SegmentVector",
                                                                                   "docstr" )
        .def( "extract_seeds", &SegmentVector::extractSeeds )
        .def( "num_seeds", &SegmentVector::numSeeds )
#if MEASURE_DURATIONS == ( 1 )
        .def_readwrite( "fExtraction", &SegmentVector::fExtraction )
        .def_readwrite( "fSorting", &SegmentVector::fSorting )
        .def_readwrite( "fLinesweep", &SegmentVector::fLinesweep )
#endif
        .def( "num_seeds_larger", &SegmentVector::numSeedsLarger );

    // tell boost python that pointers of these classes can be converted implicitly
    py::implicitly_convertible<SegmentVector, Container>( );

} // function
#endif