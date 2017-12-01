#include "getAnchors.h"


std::vector<std::shared_ptr<Container>> NlongestIntervalsAsAnchors::getInputType() const
{
    return std::vector<std::shared_ptr<Container>>{
            std::shared_ptr<Container>(new SegmentTree()),
            std::shared_ptr<Container>(new BWACompatiblePackedNucleotideSequencesCollection()),
            std::shared_ptr<Container>(new FM_Index())
        };
}//function

std::shared_ptr<Container> NlongestIntervalsAsAnchors::getOutputType() const
{
    return std::shared_ptr<Container>(new Seeds());
}//function


std::shared_ptr<Container> NlongestIntervalsAsAnchors::execute(
        std::vector<std::shared_ptr<Container>> vpInput
    )
{

    std::shared_ptr<SegmentTree> pCastedInput = std::static_pointer_cast<SegmentTree>(vpInput[0]);
	std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pRefSeq = 
		std::static_pointer_cast<BWACompatiblePackedNucleotideSequencesCollection>(vpInput[1]);
	std::shared_ptr<FM_Index> pxFM_index = std::static_pointer_cast<FM_Index>(vpInput[2]);

    std::vector<Seed> aSeeds;
    /*
    *   get the n longest intervals
    */
    pCastedInput->forEach(
        [&](std::shared_ptr<SegmentTreeInterval> pxNode)
        {
            pxNode->forEachSeed(
                pxFM_index, uiMaxHitsPerInterval, true,
                [&](Seed xS)
                {
                    aSeeds.push_back(xS);
                }//lambda
            );
        }//lambda
    );//forEach

    std::sort(
        aSeeds.begin(), aSeeds.end(),
        []
        (const Seed& a, const Seed& b)
        {
            return a.size() > b.size();
        }//lambda
    );//sort function call
    assert(aSeeds.size() <= 1 || aSeeds.front().size() >= aSeeds.back().size());

    std::shared_ptr<Seeds> pRet(new Seeds());

    for(unsigned int i = 0; i < uiN && i < aSeeds.size(); i++)
    {
        pRet->push_back(aSeeds[i]);
    }//for

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
        boost::python::init<unsigned int, unsigned int>(
            "arg1: self\n"
            "arg2: number of intervals to extract as anchors\n"
        )
    )
        .def_readwrite("uiN", &NlongestIntervalsAsAnchors::uiN)
        .def_readwrite("uiMaxHitsPerInterval", &NlongestIntervalsAsAnchors::uiMaxHitsPerInterval)
    ;

	boost::python::implicitly_convertible< 
		std::shared_ptr<NlongestIntervalsAsAnchors>,
		std::shared_ptr<CppModule> 
	>();
}//function