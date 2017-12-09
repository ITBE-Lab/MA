#include "getAnchors.h"


ContainerVector NlongestIntervalsAsAnchors::getInputType() const
{
    return ContainerVector{
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
        ContainerVector vpInput
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

    /*
     * sort them
     */
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

    // only save the uiN longest
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
        boost::python::init<unsigned int, unsigned int>()
    )
        .def_readwrite("uiN", &NlongestIntervalsAsAnchors::uiN)
        .def_readwrite("uiMaxHitsPerInterval", &NlongestIntervalsAsAnchors::uiMaxHitsPerInterval)
    ;

	boost::python::implicitly_convertible< 
		std::shared_ptr<NlongestIntervalsAsAnchors>,
		std::shared_ptr<CppModule> 
	>();
}//function