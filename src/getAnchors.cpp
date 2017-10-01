#include "getAnchors.h"


std::vector<ContainerType> NlongestIntervalsAsAnchors::getInputType()
{
    return std::vector<ContainerType>{ContainerType::segmentList};
}//function

ContainerType NlongestIntervalsAsAnchors::getOutputType()
{
    return ContainerType::segmentList;
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
            while (pxIterator.isListElement() && pxIterator->size() > pxNode->size())
                ++pxIterator;
            pRet->insertBefore(pxNode, pxIterator);
            if (pRet->length() > uiN)
                pRet->removeNode(pRet->end());
        }//lambda
    );//forEach

    return pRet;

}//function

void exportGetAnchors()
{
     //export the segmentation class
    boost::python::class_<
        NlongestIntervalsAsAnchors, 
        boost::python::bases<Module>
    >(
        "NlongestIntervalsAsAnchors", 
        "uiN: number of intervals to extract as anchors\n"
        "\n"
        "Picks a set of anchors for the strips of consideration.\n"
        "\n"
        "Execution:\n"
        "   Expects seg_list\n"
        "       seg_list: the list of segments to pick the anchors from\n"
        "   returns segList.\n"
        "       seg_list: the anchors\n",
        boost::python::init<boost::python::optional<unsigned int>>(
            "arg1: self\n"
            "arg2: number of intervals to extract as anchors\n"
        )
    )
        .add_property("uiN", &NlongestIntervalsAsAnchors::uiN, &NlongestIntervalsAsAnchors::uiN)
    ;

}//function