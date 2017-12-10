#include "segmentList.h"

void SegmentListInterval::push_back(Segment xSegment)
{
	lxSegment.push_back(xSegment);
	DEBUG(
		std::cout << "found segment: " << xSegment.saInterval().start() << " " << xSegment.saInterval().end() << std::endl;
	)
}//function

std::ostream& operator<<(std::ostream& xOs, const SegmentList& xTree)
{
	xTree.print(xOs);
	return xOs;
}//function


std::ostream& operator<<(std::ostream& xOs, const SegmentListInterval &rxNode)
{
	rxNode.print(xOs);
	return xOs;
}//function

void exportIntervalTree()
{
	 //export the SegmentListInterval class
	boost::python::class_<
			SegmentListInterval, 
			boost::python::bases<Container>,
			std::shared_ptr<SegmentListInterval>
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
				&SegmentListInterval::start_boost1,
				"arg1: self\n"
				"returns: start index of the interval\n"
			)
		.def(
				"start", 
				&SegmentListInterval::start_boost2,
				"arg1: self\n"
				"arg1: new start of the interval\n"
				"returns: nil\n"
			)
		.def(
				"end", 
				&SegmentListInterval::end_boost1,
				"arg1: self\n"
				"returns: end index of the interval\n"
			)
		.def(
				"end", 
				&SegmentListInterval::end_boost2,
				"arg1: self\n"
				"arg1: new end of the interval\n"
				"returns: nil\n"
			)
		.def(
				"push_back_bwt", 
				&SegmentListInterval::push_back,
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
				&SegmentListInterval::getSeeds,
				"arg1: self\n"
				"arg2: the fm_index.\n"
				"returns: all seeds within the segment.\n"
			)
		.def(
				"size",
				&SegmentListInterval::size_boost1,
				"arg1: self\n"
				"returns: lenght of the segment\n"
			)
		.def(
				"size",
				&SegmentListInterval::size_boost2,
				"arg1: self\n"
				"arg1: new size of the interval\n"
				"returns: nil\n"
			)
		.def(
				"set", 
				&SegmentListInterval::set,
				"arg1: self\n"
				"arg2: new start index of the interval\n"
				"arg3: new end index of the interval\n"
			);

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< 
			std::shared_ptr<SegmentListInterval>,
			std::shared_ptr<Container> 
		>();  
	
	 //export the SegmentList class
	boost::python::class_<
			SegmentList, 
			boost::python::bases<Container>,
			std::shared_ptr<SegmentList>
		>(
			"SegmentList"
		)
			.def(boost::python::init<const nucSeqIndex>())
			.def(
					"__iter__", 
					&SegmentList::begin_boost
				)
			.def(
					"get_seeds", 
					&SegmentList::getSeeds,
					boost::python::with_custodian_and_ward_postcall<1,0>(),
					"arg1: self\n"
					"arg2: the fm_index.\n"
					"returns: all seeds.\n"
				)
			.def(
					"num_seeds", 
					&SegmentList::numSeeds,
					"arg1: self\n"
					"returns: number of seeds.\n"
				)
	;

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< 
			std::shared_ptr<SegmentList>,
			std::shared_ptr<Container> 
		>(); 
	
	 //export the SegmentList iterator class
	boost::python::class_<
			SegmentList::PythonIterator
		>(
			"SegmentListIterator",
			boost::python::no_init
		)
		.def(
				"__next__", 
				&SegmentList::PythonIterator::next
			)
	;

}//function