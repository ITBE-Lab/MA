#include "module/getAnchors.h"
using namespace libLAuS;

ContainerVector GetAnchors::getInputType() const
{
    return ContainerVector{
            std::shared_ptr<Container>(new SegmentVector()),
            std::shared_ptr<Container>(new Pack()),
            std::shared_ptr<Container>(new FMIndex())
        };
}//function

std::shared_ptr<Container> GetAnchors::getOutputType() const
{
    return std::shared_ptr<Container>(new Seeds());
}//function


std::shared_ptr<Container> GetAnchors::execute(
        ContainerVector vpInput
    )
{

    std::shared_ptr<SegmentVector> pCastedInput = std::static_pointer_cast<SegmentVector>(vpInput[0]);
	std::shared_ptr<Pack> pRefSeq = 
		std::static_pointer_cast<Pack>(vpInput[1]);
	std::shared_ptr<FMIndex> pxFM_index = std::static_pointer_cast<FMIndex>(vpInput[2]);

    std::vector<Seed> aSeeds;
    /*
    *   get the n longest intervals
    */
    pCastedInput->forEachSeed(
        pxFM_index, uiMaxHitsPerInterval, true,
        [&](Seed xS)
        {
            aSeeds.push_back(xS);
        }//lambda
    );

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
        GetAnchors, 
        boost::python::bases<Module>,
		std::shared_ptr<GetAnchors>
    >(
        "GetAnchors",
        boost::python::init<unsigned int, unsigned int>()
    )
        .def_readwrite("uiN", &GetAnchors::uiN)
        .def_readwrite("uiMaxHitsPerInterval", &GetAnchors::uiMaxHitsPerInterval)
    ;

	boost::python::implicitly_convertible< 
		std::shared_ptr<GetAnchors>,
		std::shared_ptr<Module> 
	>();
}//function