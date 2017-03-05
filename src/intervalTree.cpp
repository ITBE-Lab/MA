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
