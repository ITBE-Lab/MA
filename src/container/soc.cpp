/** 
 * @file segment.cpp
 * @author Markus Schmidt
 */
#include "container/soc.h"
using namespace libMA;

void exportSoC()
{
#if DEBUG_LEVEL >= 1
        boost::python::class_<
                std::pair<nucSeqIndex, nucSeqIndex>
            >("nucSeqPair", boost::python::init<nucSeqIndex, nucSeqIndex>())
            .def_readwrite("first", &std::pair<nucSeqIndex, nucSeqIndex>::first)
            .def_readwrite("second", &std::pair<nucSeqIndex, nucSeqIndex>::second)
        ;
        boost::python::class_<std::vector<std::pair<nucSeqIndex, nucSeqIndex>>
        >("nucSeqPairVector")
        .def(boost::python::vector_indexing_suite<
                std::vector<std::pair<nucSeqIndex, nucSeqIndex>>,
                /*
                *    true = noproxy this means that the content of the vector is already exposed by
                *    boost python. 
                *    if this is kept as false, Container would be exposed a second time.
                *    the two Containers would be different and not inter castable.
                */
                true
            >());
        boost::python::class_<
                SoCPriorityQueue::blub
            >("nucSeqnucSeqInt")
            .def_readwrite("first", &SoCPriorityQueue::blub::first)
            .def_readwrite("second", &SoCPriorityQueue::blub::second)
            .def_readwrite("third", &SoCPriorityQueue::blub::third)
        ;
        boost::python::class_<std::vector<SoCPriorityQueue::blub>
        >("nucSeqnucSeqIntVector")
        .def(boost::python::vector_indexing_suite<
                std::vector<SoCPriorityQueue::blub>,
                /*
                *    true = noproxy this means that the content of the vector is already exposed by
                *    boost python. 
                *    if this is kept as false, Container would be exposed a second time.
                *    the two Containers would be different and not inter castable.
                */
                true
            >());
#endif
    //export the SoCPriorityQueue class
    boost::python::class_<
            SoCPriorityQueue, 
            boost::python::bases<Container>,
            std::shared_ptr<SoCPriorityQueue>
        >("SoCPriorityQueue")
        .def(
                "empty", 
                &SoCPriorityQueue::empty
            )
        .def(
                "pop", 
                &SoCPriorityQueue::pop
            )
        .def(
                "make_heap",
                &SoCPriorityQueue::make_heap
            )
    DEBUG(
        .def_readwrite("scores", &SoCPriorityQueue::vScores)
    )// DEBUG
    ;
}//function