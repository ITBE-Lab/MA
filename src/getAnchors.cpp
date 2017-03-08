#include "getAnchors.h"


std::shared_ptr<Container> NlongestIntervalsAsAnchors::getInputType()
{
    return std::shared_ptr<DummyContainer>(new DummyContainer(ContainerType::segmentList));
}//function

std::shared_ptr<Container> NlongestIntervalsAsAnchors::getOutputType()
{
    return std::shared_ptr<DummyContainer>(new DummyContainer(ContainerType::segmentList));
}//function


std::shared_ptr<Container> NlongestIntervalsAsAnchors::execute(std::shared_ptr<Container> pInput)
{
    std::shared_ptr<SegmentTree> pRet(new SegmentTree());

    
	std::shared_ptr<SegmentTreeContainer> pCastedInput = std::static_pointer_cast<SegmentTreeContainer>(pInput);
    /*
    *   get the n longest intervals
    */
    pCastedInput->pTree->forEach(
        [&pRet, this](std::shared_ptr<SegmentTreeInterval> pxNode)
        {
            auto pxIterator = pRet->begin();
            while (pxIterator.isListElement() && pxIterator->length() > pxNode->length())
                ++pxIterator;
            pRet->insertBefore(pxNode, pxIterator);
            if (pRet->length() > uiN)
                pRet->removeNode(pRet->end());
        }//lambda
    );//forEach

    //TODO: achors need to be extracted

    return std::shared_ptr<SegmentTreeContainer>(new SegmentTreeContainer(pRet));

}//function

void exportGetAnchors()
{
     //export the segmentation class
	boost::python::class_<
        NlongestIntervalsAsAnchors, 
        boost::python::bases<Module>
    >("NlongestIntervalsAsAnchors", boost::python::init<boost::python::optional<unsigned int>>())
		.add_property("uiN", &NlongestIntervalsAsAnchors::uiN, &NlongestIntervalsAsAnchors::uiN)
    ;

}//function