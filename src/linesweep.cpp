#if 0

#include "linesweep.h"


/** each perfect match "casts a shadow" at the left and right border of the strip
 * each shadow is stored in one of these data structures.
*/
class ShadowInterval: public Interval<nucSeqIndex>{
private:
	////while swiping the interfering shadows will get stored in this list
    std::list<ShadowInterval*> lpxInterferingIntervals;
    ////the seed container this interval corresponds to
    std::shared_ptr<SeedContainer> pxSeed;
	////the total score of the interfering shadows
    unsigned int uiScoreInterfering = 0;
public:
    ShadowInterval(
            nucSeqIndex uiBegin, 
            nucSeqIndex uiEnd, 
            std::shared_ptr<SeedContainer> pxSeed
        )
            :
        Interval(uiBegin, uiEnd),
        pxSeed(pxSeed)
    {}//constructor

    void addInterferingInterval(ShadowInterval* pInterval)
    {

    }//function

    inline const Seed& operator*() const
    {
        return *pxSeed;
    }//operator
    inline const Seed& operator->() const
    {
        return *pxSeed;
    }//operator
};//class

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

std::shared_ptr<Container> LineSweep::execute(
        std::vector<std::shared_ptr<Container>> pInput
    )
{
    std::shared_ptr<NucleotideSequence> pQuerrySeq = 
    std::static_pointer_cast<NucleotideSequence>(vpInput[0]);
    std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pRefSeq =
        std::static_pointer_cast<BWACompatiblePackedNucleotideSequencesCollection>(vpInput[1]);
    std::shared_ptr<StripOfConsideration> pStrip =
        std::static_pointer_cast<StripOfConsideration>(vpInput[2]);

    //sort (increasingly) by start coordinate of the match
		std::sort(apxPerfectMatches.begin(), apxPerfectMatches.end(),
        [](const PerfectMatchContainer xA, const PerfectMatchContainer xB)
        {
            return xA.pxPerfectMatch->getPosOnQuery() < xB.pxPerfectMatch->getPosOnQuery();
        }//lambda
    );//sort function call

    return pStrip;
}//function

void exportLinesweep()
{

}//function

#endif