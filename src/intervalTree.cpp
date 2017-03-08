#include "intervalTree.h"

std::ostream& operator<<(std::ostream& xOs, const SegmentTree& xTree)
{
	xTree.print(xOs);
	return xOs;
}//function


std::ostream& operator<<(std::ostream& xOs, const SegmentTreeInterval &rxNode)
{
	rxNode.print(xOs);
	return xOs;
}//function

void exportIntervalTree()
{
	 //export the SegmentTreeInterval class
	boost::python::class_<
			SegmentContainer, 
			boost::python::bases<Container>,
			std::shared_ptr<SegmentContainer>
		>("Segment", boost::python::init<const nucSeqIndex, const nucSeqIndex>())
			.def("start", &SegmentContainer::getStartIndex)
			.def("end", &SegmentContainer::getEndIndex)
			.def("pushBackBWT", &SegmentContainer::pushBackBwtInterval)
			.def("length", &SegmentContainer::length)
			.def("setInterval", &SegmentContainer::setInterval);

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< std::shared_ptr<SegmentContainer>, std::shared_ptr<Container> >();  
    boost::python::register_ptr_to_python< std::shared_ptr<SegmentContainer> >();
	
	 //export the SegmentTree class
	boost::python::class_<
			SegmentTreeContainer, 
			boost::python::bases<Container>,
			std::shared_ptr<SegmentTreeContainer>
		>("SegmentList")
			.def(boost::python::init<const nucSeqIndex>())
			.def("begin", &SegmentTreeContainer::begin)
			.def("end", &SegmentTreeContainer::end)
	;

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< std::shared_ptr<SegmentTreeContainer>, std::shared_ptr<Container> >(); 
    boost::python::register_ptr_to_python< std::shared_ptr<SegmentTreeContainer> >();
	
	 //export the SegmentTree iterator class
	boost::python::class_<
			SegmentListIteratorContainer
		>("SegmentListIterator", boost::python::no_init)
			.def("next", &SegmentListIteratorContainer::next)
			.def("prev", &SegmentListIteratorContainer::prev)
			.def("get", &SegmentListIteratorContainer::get)
			.def("exits", &SegmentListIteratorContainer::exits)
			.def("getCopy", &SegmentListIteratorContainer::getCopy)
	;

}//function