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

    std::shared_ptr<SegmentTree> pCastedInput = std::static_pointer_cast<SegmentTree>(vpInput[0]);

    std::vector<std::shared_ptr<SegmentTreeInterval>> aIntervals;
    /*
    *   get the n longest intervals
    */
    pCastedInput->forEach(
        [&aIntervals](std::shared_ptr<SegmentTreeInterval> pxNode)
        {
            aIntervals.push_back(pxNode);
        }//lambda
    );//forEach

    std::sort(
        aIntervals.begin(), aIntervals.end(),
        []
        (const std::shared_ptr<SegmentTreeInterval> a, const std::shared_ptr<SegmentTreeInterval> b)
        {
            return a->size() > b->size();
        }//lambda
    );//sort function call
    assert(aIntervals.front()->size() >= aIntervals.back()->size());

    std::shared_ptr<SegmentTree> pRet(new SegmentTree());

    for(unsigned int i = 0; i <= uiN && i < aIntervals.size(); i++)
        pRet->push_front(aIntervals[i]);

    return pRet;

}//function

void exportGetAnchors()
{
    boost::python::class_<
        NlongestIntervalsAsAnchors, 
        boost::python::bases<CppModule>,
		std::shared_ptr<NlongestIntervalsAsAnchors>
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

	boost::python::implicitly_convertible< 
		std::shared_ptr<NlongestIntervalsAsAnchors>,
		std::shared_ptr<CppModule> 
	>();
}//function