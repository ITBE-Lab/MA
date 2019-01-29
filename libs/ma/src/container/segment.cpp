/**
 * @file segment.cpp
 * @author Markus Schmidt
 */
#include "container/segment.h"
#include "util/pybind11.h"
using namespace libMA;

#ifdef WITH_PYTHON

#ifdef BOOST_PYTHON
void exportSegment( )
{
    // export the SegmentVector class
    boost::python::class_<Segment, boost::python::bases<Container>>( "Segment" );

    // tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible<Segment, Container>( );

    // export the SegmentVector class
    boost::python::class_<SegmentVector, boost::python::bases<Container>, std::shared_ptr<SegmentVector>>(
        "SegmentVector" )
        .def( "extract_seeds", &SegmentVector::extractSeeds
              //,boost::python::with_custodian_and_ward_postcall<1,0>()
              )
        .def( "num_seeds", &SegmentVector::numSeeds )
        .def( boost::python::vector_indexing_suite<SegmentVector,
                                                   /*
                                                    *    true = noproxy: This means that the content
                                                    * of the vector is already exposed by boost
                                                    * python. If this is kept as false, Segment would
                                                    * be exposed a second time. the two Segments would
                                                    * be different and not inter castable.
                                                    */
                                                   true>( ) )
#if MEASURE_DURATIONS == ( 1 )
        .def_readwrite( "fExtraction", &SegmentVector::fExtraction )
        .def_readwrite( "fSorting", &SegmentVector::fSorting )
        .def_readwrite( "fLinesweep", &SegmentVector::fLinesweep )
#endif
        .def( "num_seeds_larger", &SegmentVector::numSeedsLarger );

    IterableConverter( ).from_python<SegmentVector>( );

    // tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible<std::shared_ptr<SegmentVector>, std::shared_ptr<Container>>( );

} // function
#else
void exportSegment( py::module& rxPyModuleId )
{
    // export the SegmentVector class
    py::class_<Segment, Container, std::shared_ptr<Segment>>( rxPyModuleId, "Segment" );

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
#endif