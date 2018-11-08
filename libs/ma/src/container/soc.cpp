/**
 * @file segment.cpp
 * @author Markus Schmidt
 */
#include "container/soc.h"
#include "util/pybind11.h"
using namespace libMA;

#ifdef WITH_PYTHON

#ifdef BOOST_PYTHON
void exportSoC( )
{
#if DEBUG_LEVEL >= 1
    boost::python::class_<std::pair<nucSeqIndex, nucSeqIndex>>( "nucSeqPair",
                                                                boost::python::init<nucSeqIndex, nucSeqIndex>( ) )
        .def_readwrite( "first", &std::pair<nucSeqIndex, nucSeqIndex>::first )
        .def_readwrite( "second", &std::pair<nucSeqIndex, nucSeqIndex>::second );
    boost::python::class_<std::vector<std::pair<nucSeqIndex, nucSeqIndex>>>( "nucSeqPairVector" )
        .def( boost::python::vector_indexing_suite<std::vector<std::pair<nucSeqIndex, nucSeqIndex>>,
                                                   /*
                                                    * true = noproxy this means that the content of
                                                    * the vector is already exposed by boost python.
                                                    * if this is kept as false, Container would be
                                                    * exposed a second time. the two Containers
                                                    * would be different and not inter castable.
                                                    */
                                                   true>( ) );
    boost::python::class_<SoCPriorityQueue::blub>( "nucSeqnucSeqInt" )
        .def_readwrite( "first", &SoCPriorityQueue::blub::first )
        .def_readwrite( "second", &SoCPriorityQueue::blub::second )
        .def_readwrite( "qCoverage", &SoCPriorityQueue::blub::qCoverage )
        .def_readwrite( "rStart", &SoCPriorityQueue::blub::rStart )
        .def_readwrite( "rEnd", &SoCPriorityQueue::blub::rEnd )
        .def_readwrite( "rStartSoC", &SoCPriorityQueue::blub::rStartSoC )
        .def_readwrite( "rEndSoC", &SoCPriorityQueue::blub::rEndSoC );
    boost::python::class_<std::vector<SoCPriorityQueue::blub>>( "nucSeqnucSeqIntVector" )
        .def( boost::python::vector_indexing_suite<std::vector<SoCPriorityQueue::blub>,
                                                   /*
                                                    * true = noproxy this means that the content of the vector is
                                                    * already exposed by boost python. if this is kept as false,
                                                    * SoCPriorityQueue would be exposed a second time. the two
                                                    * SoCPriorityQueues would be different and not inter castable.
                                                    */
                                                   true>( ) );
    boost::python::class_<std::vector<std::shared_ptr<Seeds>>>( "seedVector" )
        .def( boost::python::vector_indexing_suite<std::vector<std::shared_ptr<Seeds>>,
                                                   /*
                                                    * true = noproxy this means that the content of the vector is
                                                    * already exposed by boost python. If this is kept as false, Seeds
                                                    * would be exposed a second time. The two Seeds would be different
                                                    * and not inter castable.
                                                    */
                                                   true>( ) );
    boost::python::class_<std::vector<double>>( "doubleVector" )
        .def( boost::python::vector_indexing_suite<std::vector<double>, true>( ) );
#endif
    // export the SoCPriorityQueue class
    boost::python::class_<SoCPriorityQueue, boost::python::bases<Container>, std::shared_ptr<SoCPriorityQueue>>(
        "SoCPriorityQueue" )
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
    // tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible<std::shared_ptr<SoCPriorityQueue>, std::shared_ptr<Container>>( );
} // function
#else
void exportSoC( py::module& rxPyModuleId )
{
#if DEBUG_LEVEL >= 1
    py::class_<std::pair<nucSeqIndex, nucSeqIndex>>( rxPyModuleId, "nucSeqPair" )
        .def( py::init<nucSeqIndex, nucSeqIndex>( ) )
        .def_readwrite( "first", &std::pair<nucSeqIndex, nucSeqIndex>::first )
        .def_readwrite( "second", &std::pair<nucSeqIndex, nucSeqIndex>::second );
    py::bind_vector<std::vector<std::pair<nucSeqIndex, nucSeqIndex>>>( rxPyModuleId, "nucSeqPairVector" );

    py::class_<SoCPriorityQueue::blub>( rxPyModuleId, "nucSeqNucSeqInterval" )
        .def_readwrite( "first", &SoCPriorityQueue::blub::first )
        .def_readwrite( "second", &SoCPriorityQueue::blub::second )
        .def_readwrite( "qCoverage", &SoCPriorityQueue::blub::qCoverage )
        .def_readwrite( "rStart", &SoCPriorityQueue::blub::rStart )
        .def_readwrite( "rEnd", &SoCPriorityQueue::blub::rEnd )
        .def_readwrite( "rStartSoC", &SoCPriorityQueue::blub::rStartSoC )
        .def_readwrite( "rEndSoC", &SoCPriorityQueue::blub::rEndSoC );
    py::bind_vector<std::vector<SoCPriorityQueue::blub>>( rxPyModuleId, "nucSeqNucSeqIntervalVector" );

    py::bind_vector<std::vector<std::shared_ptr<Seeds>>>( rxPyModuleId, "seedVector" );
    py::bind_vector<std::vector<double>>( rxPyModuleId, "doubleVector" );
#endif
    // export the SoCPriorityQueue class
    py::class_<SoCPriorityQueue, Container, std::shared_ptr<SoCPriorityQueue>>( rxPyModuleId, "SoCPriorityQueue" )
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
    py::implicitly_convertible<SoCPriorityQueue, Container>( );
} // function
#endif
#endif
