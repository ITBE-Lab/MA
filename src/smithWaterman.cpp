#include "smithWaterman.h"


std::vector<ContainerType> SmithWaterman::getInputType()
{
    return std::vector<ContainerType>{
        //the sound strip of consideration
        ContainerType::segmentList,
        //the query sequence
        ContainerType:nucSeq,
        //the reference sequence
        ContainerType:packedNucSeq,
    };
}//function

std::vector<ContainerType> SmithWaterman::getOutputType()
{
    return std::vector<ContainerType>{ContainerType::alignment};
}//function


std::shared_ptr<Container> SmithWaterman::execute(
        std::vector<std::shared_ptr<Container>> vpInput
    )
{

    std::shared_ptr<SegmentTree> pCastedInput = std::static_pointer_cast<SegmentTree>(vpInput[0]);
    std::shared_ptr<NucleotideSequence> pQuery 
        = std::static_pointer_cast<NucleotideSequence>(vpInput[1]);
    std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pRef = 
        std::static_pointer_cast<BWACompatiblePackedNucleotideSequencesCollection>(vpInput[2]);

    return nullptr;

}//function

void exportSmithWaterman()
{
     //export the segmentation class
    boost::python::class_<
        SmithWaterman, 
        boost::python::bases<Module>
    >(
        "SmithWaterman", 
        "Picks a set of anchors for the strips of consideration.\n"
        "\n"
        "Execution:\n"
        "   Expects seg_list, query, ref\n"
        "       seg_list: the list of segments to pick the anchors from\n"
        "       query: the query as NucSeq\n"
        "       ref: the reference as Pack\n"
        "   returns alignment.\n"
        "       alignment: the final alignment\n"
    );

}//function