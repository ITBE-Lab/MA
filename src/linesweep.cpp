#include "linesweep.h"

std::vector<ContainerType> LineSweep::getInputType()
{
	return std::vector<ContainerType>{
			//the querry
			ContainerType::nucSeq,
			//the reference
			ContainerType::packedNucSeq,
			//the stips of consideration
			ContainerType::stripOfConsideration,
		};
}//function

std::vector<ContainerType> LineSweep::getOutputType()
{
	return std::vector<ContainerType>{ContainerType::stripOfConsideration};
}//function

/**
 * determine the start and end positions this match casts on the left border of the given bucket
 * pxMatch is the container this match is stored in.
 */
ShadowInterval LineSweep::getLeftShadow(
        nucSeqIndex uiBucketStart,
        std::list<Seed>::iterator pSeed,
        nucSeqIndex uiBucketSize,
        nucSeqIndex uiQueryLength
    ) const
{
    return ShadowInterval(
            pSeed->start() + uiBucketStart,
            pSeed->start_ref() + pSeed->size() + uiQueryLength + uiBucketSize,
            pSeed
        );
}//function

/**
 * determine the start and end positions this match casts on the right border of the given bucket
 * pxMatch is the container this match is stored in.
 */
ShadowInterval LineSweep::getRightShadow(
    nucSeqIndex iBucketStart,
    std::list<Seed>::iterator pSeed,
    nucSeqIndex iBucketSize,
    nucSeqIndex iQueryLength) const
{
    return ShadowInterval(
            pSeed->start_ref(),
            pSeed->start() + pSeed->size() + iQueryLength + iBucketSize * 2 + iBucketStart,
            pSeed
        );
}//function

void LineSweep::linesweep(std::vector<ShadowInterval>& vShadows, std::list<Seed>& rSeeds)
{
    //sort shadows (increasingly) by start coordinate of the match
    std::sort(
        vShadows.begin(),
        vShadows.end(),
        [](ShadowInterval xA, ShadowInterval xB)
        {
            return xA.start() < xB.start();
        }//lambda
    );//sort function call

    //records the interval ends
    SelfBalancingBinarySearchTree<ShadowInterval*> xItervalEnds = SelfBalancingBinarySearchTree<ShadowInterval*>();

    //this is the line sweeping part
    for(ShadowInterval& rInterval : vShadows)
    {
        //TODO: remove the stack markers
        while(true)
        {
            ShadowInterval* pFirstEnding = xItervalEnds.first();
            //check if we really need to remove the first interval in the tree
            if(pFirstEnding == nullptr)
                break;
            if(pFirstEnding->end() >= rInterval.start())
                break;

            //when reaching here we actually have to remove the intervall
            pFirstEnding->removeSeedIfNecessary(rSeeds);
            xItervalEnds.deleteFirst();
        }//while

        //work on the current interval
        ShadowInterval** ppNextShadow = xItervalEnds.insert(&rInterval);
        if(ppNextShadow != nullptr)
            (*ppNextShadow)->addInterferingInterval(&rInterval);
    }//for
}//function

std::shared_ptr<Container> LineSweep::execute(
        std::vector<std::shared_ptr<Container>> vpInput
    )
{
    std::shared_ptr<NucleotideSequence> pQuerrySeq =
        std::static_pointer_cast<NucleotideSequence>(vpInput[0]);
    std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pRefSeq =
        std::static_pointer_cast<BWACompatiblePackedNucleotideSequencesCollection>(vpInput[1]);
    std::shared_ptr<StripOfConsideration> pStrip =
        std::static_pointer_cast<StripOfConsideration>(vpInput[2]);

    std::vector<ShadowInterval> vLeftShadows = {};

    //get the left shadows
    pStrip->forAllSeeds(
        [&](std::list<Seed>::iterator pSeed)
        {
            vLeftShadows.push_back(getLeftShadow(
                    pStrip->start(), 
                    pSeed, 
                    pStrip->size(), 
                    pQuerrySeq->length()
                ));
        }//lambda
    );

    //perform the line sweep algorithm on the left shadows
    linesweep(vLeftShadows, pStrip->seeds());

    return pStrip;
}//function

void exportLinesweep()
{
    //export the LineSweepContainer class
	boost::python::class_<LineSweep, boost::python::bases<Module>>(
        "LineSweep",
        "Uses linesweeping to remove contradicting "
        "matches within one strip of consideration.\n"
        "\n"
        "Execution:\n"
        "	Expects querry, ref, strip_vec as input.\n"
        "		querry: the querry as NucleotideSequence\n"
        "		ref: the reference seqeuence as Pack\n"
        "		strip_vec: the areas that shall be evaluated as StripOfConsiderationVector\n"
        "	returns strip_vec.\n"
        "		strip_vec: the evaluated areas\n"
    );
}//function