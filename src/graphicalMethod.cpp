#include "graphicalMethod.h"


std::shared_ptr<Container> Bucketing::getInputType()
{
	std::shared_ptr<ContainerVector> pRet(new ContainerVector());
	//all segments
	pRet->vElements.push_back(std::shared_ptr<Container>(new DummyContainer(ContainerType::segmentList)));
	//the anchors
	pRet->vElements.push_back(std::shared_ptr<Container>(new DummyContainer(ContainerType::segmentList)));
	return pRet;
}//function

std::shared_ptr<Container> Bucketing::getOutputType()
{
	return std::shared_ptr<ContainerVector> pRet(new DummyContainer(ContainerType::stripOfConsiderationList));
}//function


void exportGraphicalMethod()
{
    //export the segmentation class
	boost::python::class_<LineSweepContainer, boost::python::bases<Module>>("LineSweep")
		;
}