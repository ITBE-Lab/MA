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
			SegmentTreeInterval, 
			boost::python::bases<Container>,
			std::shared_ptr<SegmentTreeInterval>
		>("Segment", boost::python::init<const nucSeqIndex, const nucSeqIndex>())
			.def("start", &SegmentTreeInterval::getStartIndex)
			.def("end", &SegmentTreeInterval::getEndIndex)
			.def("pushBackBWT", &SegmentTreeInterval::pushBackBwtInterval)
			.def("length", &SegmentTreeInterval::length)
			.def("setInterval", &SegmentTreeInterval::setInterval);

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< 
			std::shared_ptr<SegmentTreeInterval>,
			std::shared_ptr<Container> 
		>();  
    boost::python::register_ptr_to_python< std::shared_ptr<SegmentTreeInterval> >();
	
	 //export the SegmentTree class
	boost::python::class_<
			SegmentTree, 
			boost::python::bases<Container>,
			std::shared_ptr<SegmentTree>
		>("SegmentList")
			.def(boost::python::init<const nucSeqIndex>())
			.def("begin", &SegmentTree::begin)
			.def("end", &SegmentTree::end)
	;

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< 
			std::shared_ptr<SegmentTree>,
			std::shared_ptr<Container> 
		>(); 
    boost::python::register_ptr_to_python< std::shared_ptr<SegmentTree> >();
	
	 //export the SegmentTree iterator class
	boost::python::class_<
			SegmentTree::Iterator
		>("SegmentListIterator", boost::python::no_init)
			.def("next", &SegmentTree::Iterator::operator++)
			.def("prev", &SegmentTree::Iterator::operator--)
			.def("get", &SegmentTree::Iterator::operator*)
			.def("exits", &SegmentTree::Iterator::isListElement)
			.def("getCopy", &SegmentTree::Iterator::getCopy)
	;

}//function