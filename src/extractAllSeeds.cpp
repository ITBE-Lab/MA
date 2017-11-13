#include "extractAllSeeds.h"



std::vector<ContainerType> ExtractAllSeeds::getInputType()
{
	return std::vector<ContainerType>
	{
		//all segments
		ContainerType::segmentList,
		//the forward fm_index
		ContainerType::fM_index
	};
}//function

ContainerType ExtractAllSeeds::getOutputType()
{
	return ContainerType::seeds;
}//function



std::shared_ptr<Container> ExtractAllSeeds::execute(
		std::vector<std::shared_ptr<Container>> vpInput
	)
{
	std::shared_ptr<SegmentTree> pSegments = std::static_pointer_cast<SegmentTree>(vpInput[0]);
	std::shared_ptr<FM_Index> pFM_index = std::static_pointer_cast<FM_Index>(vpInput[1]);


	return pSegments->getSeeds(pFM_index, maxNum);
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