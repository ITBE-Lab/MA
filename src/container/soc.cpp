/** 
 * @file segment.cpp
 * @author Markus Schmidt
 */
#include "container/soc.h"
using namespace libMA;

void exportSoC()
{
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
    ;
}//function