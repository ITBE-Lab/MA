#include "module/getAnchors.h"
using namespace libMABS;

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
        std::shared_ptr<ContainerVector> vpInput
    )
{
    std::shared_ptr<SegmentVector> pCastedInput =
        std::static_pointer_cast<SegmentVector>((*vpInput)[0]);
    std::shared_ptr<Pack> pRefSeq = 
        std::static_pointer_cast<Pack>((*vpInput)[1]);
    std::shared_ptr<FMIndex> pxFM_index = std::static_pointer_cast<FMIndex>((*vpInput)[2]);

    std::shared_ptr<Seeds> pRet = std::shared_ptr<Seeds>(new Seeds());
    std::set<nucSeqIndex> xTakenSpots;
    /*
    * sort the intervals on the query by their length
    */
    if(!bAlreadySorted)
        std::sort(
            begin(), end(),
            []
            (const std::shared_ptr<Segment> a, std::shared_ptr<Segment> b)
            {
                return a->size() > b->size();
            }//lambda
        );//sort function call
    assert(size() <= 1 || front()->size() >= back()->size());
    //@todo maybe we don't want exclude every anchor in the entire strip
    nucSeqIndex uiStripSize = StripOfConsideration::getStripSize(pQuerySeq->length());
    pCastedInput->forEachSeed(
        pxFMIndex,
        uiMAxAmbiguity,
        bSkip,
        [&pRet, &uiN]
        (Seed s)
        {   
            //+ uiStripSize to avoid underflows
            nucSeqIndex uiRootPos = StripOfConsideration::getPositionForBucketing(pQuerySeq->length(), s) + uiStripSize;

            auto xIterator = xTakenSpots.lower_bound(uiRootPos - uiStripSize);
            if(xIterator != xTakenSpots.end() && *xIterator <= uiRootPos + uiStripSize)
                //continue cause this anchor will already be collected in a previous strip
                return true;

            //@todo check uppter bound and add uiRootPos to xTakenSpots

            //@todo filter here
            pRet->push_back(s);
            if(pRet->size() >= uiN)
                return false;
            return true;
        }//lambda
    );//for each
    assert(pRet->size() <= 1 || pRet->front().size() >= pRet->back().size());
    return pRet;

    return pCastedInput->extractLargestSeeds(pxFM_index, uiN, uiMaxAmbiguity);
}//function

void exportGetAnchors()
{
    boost::python::class_<
        GetAnchors, 
        boost::python::bases<Module>,
        std::shared_ptr<GetAnchors>
    >("GetAnchors")
        .def(boost::python::init<unsigned int, unsigned int>())
        .def_readwrite("n", &GetAnchors::uiN)
        .def_readwrite("max_ambiguity", &GetAnchors::uiMaxAmbiguity)
    ;

    boost::python::implicitly_convertible< 
        std::shared_ptr<GetAnchors>,
        std::shared_ptr<Module> 
    >();
}//function