#include "intervalTree.h"

void SegmentTreeInterval::push_back(SaSegment xSaSegment)
{
	lxSaSegment.push_back(xSaSegment);
	DEBUG(
		std::cout << "found segment: " << xSaSegment.saInterval().start() << " " << xSaSegment.saInterval().end() << std::endl;
	)
}//function

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
		>(
			"Segment",
			"A interval on the querry.\n",
			boost::python::init<const nucSeqIndex, const nucSeqIndex>(
				"arg1: self\n"
				"arg2: start index of the segment\n"
				"arg2: size of the segment\n"
			)
		)
		.def(
				"start", 
				&SegmentTreeInterval::start_boost1,
				"arg1: self\n"
				"returns: start index of the interval\n"
			)
		.def(
				"start", 
				&SegmentTreeInterval::start_boost2,
				"arg1: self\n"
				"arg1: new start of the interval\n"
				"returns: nil\n"
			)
		.def(
				"end", 
				&SegmentTreeInterval::end_boost1,
				"arg1: self\n"
				"returns: end index of the interval\n"
			)
		.def(
				"end", 
				&SegmentTreeInterval::end_boost2,
				"arg1: self\n"
				"arg1: new end of the interval\n"
				"returns: nil\n"
			)
		.def(
				"push_back_bwt", 
				&SegmentTreeInterval::push_back,
				"arg1: self\n"
				"arg2: start index of the BWT Interval\n"
				"arg3: length of the BWT Interval\n"
				"arg4: start index on the query\n"
				"arg5: end index on the query\n"
				"arg6: weather the BWT Interval is on the forwars or reversed FM Index\n"
				"arg7: shall the hit be used as anchor for strips of consideration\n"
				"returns: nil\n"
				"\n"
				"Record a bwt interval that matches this segment (or a part of this segment).\n"
				"use uiStartOfIntervalOnQuery and uiEndOfIntervalOnQuery to record the actual "
				"size of the match\n"
			)
		.def(
				"get_seeds", 
				&SegmentTreeInterval::getSeeds,
				"arg1: self\n"
				"arg2: the fm_index.\n"
				"returns: all seeds within the segment.\n"
			)
		.def(
				"size",
				&SegmentTreeInterval::size_boost1,
				"arg1: self\n"
				"returns: lenght of the segment\n"
			)
		.def(
				"size",
				&SegmentTreeInterval::size_boost2,
				"arg1: self\n"
				"arg1: new size of the interval\n"
				"returns: nil\n"
			)
		.def(
				"set", 
				&SegmentTreeInterval::set,
				"arg1: self\n"
				"arg2: new start index of the interval\n"
				"arg3: new end index of the interval\n"
			);

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< 
			std::shared_ptr<SegmentTreeInterval>,
			std::shared_ptr<Container> 
		>();  
	
	 //export the SegmentTree class
	boost::python::class_<
			SegmentTree, 
			boost::python::bases<Container>,
			std::shared_ptr<SegmentTree>
		>(
			"SegmentList"
		)
			.def(boost::python::init<const nucSeqIndex>())
			.def(
					"__iter__", 
					&SegmentTree::begin_boost
				)
			.def(
					"get_seeds", 
					&SegmentTree::getSeeds,
					boost::python::with_custodian_and_ward_postcall<1,0>(),
					"arg1: self\n"
					"arg2: the fm_index.\n"
					"returns: all seeds.\n"
				)
			.def(
					"num_seeds", 
					&SegmentTree::numSeeds,
					"arg1: self\n"
					"returns: number of seeds.\n"
				)
	;

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< 
			std::shared_ptr<SegmentTree>,
			std::shared_ptr<Container> 
		>(); 
	
	 //export the SegmentTree iterator class
	boost::python::class_<
			SegmentTree::PythonIterator
		>(
			"SegmentListIterator",
			boost::python::no_init
		)
		.def(
				"__next__", 
				&SegmentTree::PythonIterator::next
			)
	;

}//function