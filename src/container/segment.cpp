/** 
 * @file segment.cpp
 * @author Markus Schmidt
 */
#include "container/segment.h"
using namespace libMA;

#ifdef WITH_PYTHON
void exportIntervalTree()
{
    //export the SegmentVector class
    boost::python::class_<
            Segment, 
            boost::python::bases<Container>,
            std::shared_ptr<Segment>
        >("Segment");

    //tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible<
            std::shared_ptr<Segment>,
            std::shared_ptr<Container> 
        >();

    //export the SegmentVector class
    boost::python::class_<
            SegmentVector, 
            boost::python::bases<Container>,
            std::shared_ptr<SegmentVector>
        >("SegmentVector")
            .def(
                    "extract_seeds", 
                    &SegmentVector::extractSeeds
                    //,boost::python::with_custodian_and_ward_postcall<1,0>()
                )
            .def(
                    "num_seeds", 
                    &SegmentVector::numSeeds
                )
            .def(boost::python::vector_indexing_suite
                    <
                        SegmentVector,
                        /*
                        *    true = noproxy: This means that the content of
                        *    the vector is already exposed by boost python.
                        *    If this is kept as false, Segment would be exposed a second time.
                        *    the two Segments would be different and not inter castable.
                        */
                        true
                    >()
                )
#if MEASURE_DURATIONS == ( 1 )
            .def_readwrite("fExtraction", &SegmentVector::fExtraction)
            .def_readwrite("fSorting", &SegmentVector::fSorting)
            .def_readwrite("fLinesweep", &SegmentVector::fLinesweep)
#endif
    ;

    IterableConverter()
        .from_python<SegmentVector>();

    //tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible<
            std::shared_ptr<SegmentVector>,
            std::shared_ptr<Container> 
        >();

}//function
#endif