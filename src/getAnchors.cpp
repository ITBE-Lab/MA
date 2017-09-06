#include "getAnchors.h"


std::vector<ContainerType> NlongestIntervalsAsAnchors::getInputType()
{
    return std::vector<ContainerType>{ContainerType::segmentList};
}//function

std::vector<ContainerType> NlongestIntervalsAsAnchors::getOutputType()
{
    return std::vector<ContainerType>{ContainerType::segmentList};
}//function


std::shared_ptr<Container> NlongestIntervalsAsAnchors::execute(
        std::vector<std::shared_ptr<Container>> vpInput
    )
{
    std::shared_ptr<SegmentTree> pRet(new SegmentTree());

    std::shared_ptr<SegmentTree> pCastedInput = std::static_pointer_cast<SegmentTree>(vpInput[0]);
    /*
    *   get the n longest intervals
    */
    pCastedInput->forEach(
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

    return pRet;

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