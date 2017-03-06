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
			boost::noncopyable, 
			boost::python::bases<Container>, 
			std::shared_ptr<SegmentTreeInterval>
		>("Segment", boost::python::init<const nucSeqIndex, const nucSeqIndex>())
			.def("start", &SegmentTreeInterval::getStartIndex)
			.def("end", &SegmentTreeInterval::getEndIndex)
			.def("pushBackBWT", &SegmentTreeInterval::pushBackBwtInterval)
			.def("length", &SegmentTreeInterval::length)
			.def("setInterval", &SegmentTreeInterval::setInterval);

	//tell boost python that it's possible to convert shared pointers with these classes
    boost::python::implicitly_convertible<std::shared_ptr<SegmentTreeInterval>,std::shared_ptr<Container>>();

	
	 //export the SegmentTree class
	boost::python::class_<
			SegmentTree, 
			boost::noncopyable, 
			boost::python::bases<Container>, 
			std::shared_ptr<SegmentTree>
		>("SegmentList")
			.def(boost::python::init<const nucSeqIndex>())
			.def("getAnchors", &SegmentTree::getTheNLongestIntervals);

}//function