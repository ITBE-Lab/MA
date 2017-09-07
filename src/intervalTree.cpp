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
		>(
			"Segment",
			"A interval on the querry.\n",
			boost::python::init<const nucSeqIndex, const nucSeqIndex>()
		)
			.def(
					"start", 
					&SegmentTreeInterval::getStartIndex,
					"	arg1: self\n"
					"	returns: start index of the interval\n"
				)
			.def(
					"end", 
					&SegmentTreeInterval::getEndIndex,
					"	arg1: self\n"
					"	returns: end index of the interval\n"
				)
			.def(
					"push_back_bwt", 
					&SegmentTreeInterval::pushBackBwtInterval,
					"	arg1: self\n"
					"	arg2: start index of the BWT Interval\n"
					"	arg3: length of the BWT Interval\n"
					"	arg4: start index on the query\n"
					"	arg5: end index on the query\n"
					"	arg6: weather the BWT Interval is on the forwars or reversed FM Index\n"
					"	arg7: shall the hit be used as anchor for strips of consideration\n"
					"	returns: nil\n"
					"\n"
					"Record a bwt interval that matches this segment (or a part of this segment).\n"
					"use uiStartOfIntervalOnQuery and uiEndOfIntervalOnQuery to record the actual "
					"size of the match\n"
				)
			.def(
					"length", 
					&SegmentTreeInterval::length,
					"	arg1: self\n"
					"	returns: lenght of the segment\n"
				)
			.def(
					"set_interval", 
					&SegmentTreeInterval::setInterval,
					"	arg1: self\n"
					"	arg2: new start index of the interval\n"
					"	arg3: new end index of the interval\n"
				);

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
		>(
			"SegmentList",
			"A doubly linked list holding Sements.\n"
		)
			.def(boost::python::init<const nucSeqIndex>())
			.def(
					"begin", 
					&SegmentTree::begin,
					"	arg1: self\n"
					"	returns: an itterator pointing to the first element of the list.\n"
				)
			.def(
					"end", 
					&SegmentTree::end,
					"	arg1: self\n"
					"	returns: an itterator pointing to the last element of the list.\n"
				)
			.def(
					"insert_before", 
					&SegmentTree::insertBefore_boost,
					"	arg1: self\n"
					"	arg2: the segment to insert.\n"
					"	arg3 the current iterator.\n"
					"	returns: a new iterator.\n"
					"\n"
					"Insert a segment before the position of the given iterator.\n"
					"Invalidates the iterator.\n"
					"Replace it with the returned one.\n"
				)
			.def(
					"insert_after", 
					&SegmentTree::insertAfter_boost,
					"	arg1: self\n"
					"	arg2: the segment to insert.\n"
					"	arg3: the current iterator.\n"
					"	returns: a new iterator.\n"
					"\n"
					"Insert a segment after the position of the given iterator.\n"
					"Invalidates the iterator.\n"
					"Replace it with the returned one.\n"
				)
			.def(
					"push_back", 
					&SegmentTree::push_back,
					"	arg1: self\n"
					"	arg2: the segment to insert.\n"
					"	returns: a new iterator.\n"
					"\n"
					"Insert a segment at the end of the list.\n"
				)
			.def(
					"push_front", 
					&SegmentTree::push_front,
					"	arg1: self\n"
					"	arg2: the segment to insert.\n"
					"	returns: a new iterator.\n"
					"\n"
					"Insert a segment at the front of the list.\n"
				)
			.def(
					"remove_node", 
					&SegmentTree::removeNode_boost,
					"	arg1: self\n"
					"	arg2: the iterator pointing to the segment that shall be removed.\n"
					"	returns: nil.\n"
					"\n"
					"Invalidates the iterator.\n"
				)
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
		>(
			"SegmentListIterator",
			"Itterator for the doubely linked SegmentList.\n"
			"The list has a dommy element at bot ends.\n"
			"If the itterator sits ontop a dummy element exists() returns false.\n",
			boost::python::no_init
		)
			.def(
					"next", 
					&SegmentTree::Iterator::operator++,
					"	arg1: self\n"
					"	returns: self\n"
					"\n"
					"Move the itterator to the next element of the list.\n"
				)
			.def(
					"prev", 
					&SegmentTree::Iterator::operator--,
					"	arg1: self\n"
					"	returns: self\n"
					"\n"
					"Move the itterator to the previous element of the list.\n"
				)
			.def(
					"get", 
					&SegmentTree::Iterator::operator*,
					"	arg1: self\n"
					"	returns: the current Segment\n"
					"\n"
					"Get the Segment the itterator is currently sitting on.\n"
				)
			.def(
					"exists", 
					&SegmentTree::Iterator::isListElement,
					"	arg1: self\n"
					"	returns: weather the itterator sits on a actual element of the List.\n"
				)
			.def(
					"getCopy", 
					&SegmentTree::Iterator::getCopy,
					"	arg1: self\n"
					"	returns: copy of the itterator.\n"
				)
	;

}//function