/**
 * @file segment.cpp
 * @author Markus Schmidt
 */
#include "ma/container/soc.h"
#include "ms/util/pybind11.h"
using namespace libMA;

#ifdef WITH_PYTHON

void exportSoC( libMS::SubmoduleOrganizer& xOrganizer )
{
#if DEBUG_LEVEL >= 1
    py::bind_vector<std::vector<std::pair<nucSeqIndex, nucSeqIndex>>>( xOrganizer._util(), "nucSeqPairVector", "docstr" );

    py::class_<SoCPriorityQueue::blub>( xOrganizer._util(), "nucSeqNucSeqInterval" )
        .def_readwrite( "first", &SoCPriorityQueue::blub::first )
        .def_readwrite( "second", &SoCPriorityQueue::blub::second )
        .def_readwrite( "qCoverage", &SoCPriorityQueue::blub::qCoverage )
        .def_readwrite( "rStart", &SoCPriorityQueue::blub::rStart )
        .def_readwrite( "rEnd", &SoCPriorityQueue::blub::rEnd )
        .def_readwrite( "rStartSoC", &SoCPriorityQueue::blub::rStartSoC )
        .def_readwrite( "rEndSoC", &SoCPriorityQueue::blub::rEndSoC );
    py::bind_vector<std::vector<SoCPriorityQueue::blub>>( xOrganizer._util(), "nucSeqNucSeqIntervalVector", "docstr" );
#endif
    py::bind_vector<std::vector<std::shared_ptr<Seeds>>>( xOrganizer.util(), "seedVector", "docstr" );
    // export the SoCPriorityQueue class
    py::class_<SoCPriorityQueue, libMS::Container, std::shared_ptr<SoCPriorityQueue>>( xOrganizer.container(),
                                                                                       "SoCPriorityQueue" )
        .def( py::init<>( ) )
        .def( "empty", &SoCPriorityQueue::empty )
        .def( "pop", &SoCPriorityQueue::pop )
        .def( "make_heap", &SoCPriorityQueue::make_heap )
        .def( "__len__", &SoCPriorityQueue::size )
            DEBUG(.def_readwrite( "scores", &SoCPriorityQueue::vScores )
                      .def_readwrite( "extract", &SoCPriorityQueue::vExtractOrder )
                      .def_readwrite( "vSoCs", &SoCPriorityQueue::vSoCs )
                      .def_readwrite( "vHarmSoCs", &SoCPriorityQueue::vHarmSoCs )
                      .def_readwrite( "vSlopes", &SoCPriorityQueue::vSlopes )
                      .def_readwrite( "vIntercepts", &SoCPriorityQueue::vIntercepts )
                      .def_readwrite( "vIngroup", &SoCPriorityQueue::vIngroup ) ) // DEBUG
        ;
    // tell python that pointers of these classes can be converted implicitly
    py::implicitly_convertible<SoCPriorityQueue, libMS::Container>( );
} // function
#endif
