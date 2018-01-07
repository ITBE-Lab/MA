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