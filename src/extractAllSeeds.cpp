#include "extractAllSeeds.h"



ContainerVector ExtractAllSeeds::getInputType() const
{
	return ContainerVector
	{
		//all segments
		std::shared_ptr<Container>(new SegmentVector()),
		//the forward fm_index
		std::shared_ptr<Container>(new FMIndex())
	};
}//function

std::shared_ptr<Container> ExtractAllSeeds::getOutputType() const
{
	return std::shared_ptr<Container>(new Seeds());
}//function



std::shared_ptr<Container> ExtractAllSeeds::execute(
		ContainerVector vpInput
	)
{
	std::shared_ptr<SegmentVector> pSegments = std::static_pointer_cast<SegmentVector>(vpInput[0]);
	std::shared_ptr<FMIndex> pFM_index = std::static_pointer_cast<FMIndex>(vpInput[1]);


	return pSegments->extractSeeds(pFM_index, maxNum);
}//function

void exportExtractAllSeeds()
{
    //export the ExtractAllSeeds class
	boost::python::class_<
			ExtractAllSeeds, 
			boost::python::bases<CppModule>, 
            std::shared_ptr<ExtractAllSeeds>
		>(
        "ExtractAllSeeds"
    )
        .def_readwrite("max_hits", &ExtractAllSeeds::maxNum);

	boost::python::implicitly_convertible< 
		std::shared_ptr<ExtractAllSeeds>,
		std::shared_ptr<CppModule> 
	>();
}//function