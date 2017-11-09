#include "linesweep.h"

std::vector<ContainerType> LineSweep::getInputType()
{
	return std::vector<ContainerType>{
			ContainerType::seeds
		};
}//function

ContainerType LineSweep::getOutputType()
{
	return ContainerType::seeds;
}//function

/**
 * determine the start and end positions this match casts on the left border of the given bucket
 * pxMatch is the container this match is stored in.
 */
ShadowInterval LineSweep::getLeftShadow(std::list<Seed>::iterator pSeed) const
{
    return ShadowInterval(
            pSeed->start(),
            pSeed->end_ref() - (int64_t)pSeed->start(),
            pSeed
        );
}//function

/**
 * determine the start and end positions this match casts on the right border of the given bucket
 * pxMatch is the container this match is stored in.
 */
ShadowInterval LineSweep::getRightShadow(std::list<Seed>::iterator pSeed) const
{
    return ShadowInterval(
            pSeed->start_ref(),
            pSeed->end() - (int64_t)pSeed->start_ref(),
            pSeed
        );
}//function

void LineSweep::linesweep(
        std::vector<ShadowInterval>& vShadows, 
        std::shared_ptr<Seeds> pSeeds
    )
{
    //sort shadows (increasingly) by start coordinate of the match
    std::sort(
            vShadows.begin(),
            vShadows.end(),
            [](ShadowInterval xA, ShadowInterval xB)
            {
                /*
                * sort by the interval starts
                * if two intervals start at the same point the larger one shall be treated first
                */
                if(xA.start() == xB.start())
                    return xA.end() > xB.end();
                return xA.start() < xB.start();
            }//lambda
        );//sort function call

    //records the interval ends
    SearchTree<ShadowIntervalPtr> xItervalEnds;

    //this is the line sweeping part
    for(ShadowInterval& rInterval : vShadows)
    {
        //remove the stack markers
        while(!xItervalEnds.isEmpty())
        {

            ShadowIntervalPtr pFirstEnding = xItervalEnds.first();
            //check if we really need to remove the first interval in the tree
            if(pFirstEnding->end() > rInterval.start())
                break;

            DEBUG(
                std::cout << "Current Sweep position: " << pFirstEnding->end() << std::endl;
                std::cout << "\tend of interval " << pFirstEnding->start() << ", "
                    << pFirstEnding->end() << std::endl;
            std::cout << "\t(start_ref, end_ref; start_query, end_query) " 
                << pFirstEnding->pSeed->start_ref() << ", " << pFirstEnding->pSeed->end_ref() << "; "
                << pFirstEnding->pSeed->start() << ", " << pFirstEnding->pSeed->end() << std::endl;
            )

            //when reaching here we actually have to remove the intervall
            pFirstEnding->removeSeedIfNecessary(pSeeds);
            xItervalEnds.deleteFirst();
        }//while

        DEBUG(
            std::cout << "Current Sweep position: " << rInterval.start() << std::endl;
            std::cout << "\tat start of interval " << rInterval.start() <<
                ", " << rInterval.end() << std::endl;
            std::cout << "\t(start_ref, end_ref; start_query, end_query) " 
                << rInterval.pSeed->start_ref() << ", " << rInterval.pSeed->end_ref() << "; "
                << rInterval.pSeed->start() << ", " << rInterval.pSeed->end() << std::endl;
        )

        SearchTree<ShadowIntervalPtr>::Iterator pNextShadow = xItervalEnds.insert(
                ShadowIntervalPtr(&rInterval)
            );
        ++pNextShadow;
        if(pNextShadow.exists())
        {
            pNextShadow->addInterferingInterval(rInterval);
            SearchTree<ShadowIntervalPtr>::Iterator pFollowingShadows = pNextShadow;
            ++pFollowingShadows;
            while(
                    pFollowingShadows.exists() && 
                    //&* converts from iterator to ptr
                    (*pFollowingShadows).p != pNextShadow->getIInterferWith()
                )
            {
                pFollowingShadows->add2ndOrderInterferingInterval(rInterval);
                ++pFollowingShadows;
            }//while
        }//if
        DEBUG(
            else
                std::cout << "\tno interference" << std::endl;
        )
    }//for
    
    //remove the stack markers
    while(!xItervalEnds.isEmpty())
    {

        ShadowIntervalPtr pFirstEnding = xItervalEnds.first();

        DEBUG(
            std::cout << "Current Sweep position: " << pFirstEnding->end() << std::endl;
            std::cout << "\tat end of interval (cleanup) " << pFirstEnding->start() << ", " 
                << pFirstEnding->end() << std::endl;
            std::cout << "\t(start_ref, end_ref; start_query, end_query) " 
                << pFirstEnding->pSeed->start_ref() << ", " << pFirstEnding->pSeed->end_ref() << "; "
                << pFirstEnding->pSeed->start() << ", " << pFirstEnding->pSeed->end() << std::endl;
        )

        //when reaching here we actually have to remove the intervall
        pFirstEnding->removeSeedIfNecessary(pSeeds);
        xItervalEnds.deleteFirst();
    }//while

}//function

std::shared_ptr<Container> LineSweep::execute(
        std::vector<std::shared_ptr<Container>> vpInput
    )
{
    std::shared_ptr<Seeds> pSeeds = std::shared_ptr<Seeds>(new Seeds(
        std::static_pointer_cast<Seeds>(vpInput[0])));

    std::vector<ShadowInterval> vShadows = {};

    //get the left shadows
    for(std::list<Seed>::iterator pSeed = pSeeds->begin(); pSeed != pSeeds->end(); pSeed++)
        vShadows.push_back(getLeftShadow(pSeed));

    //perform the line sweep algorithm on the left shadows
    linesweep(vShadows, pSeeds);
    
    vShadows.clear();

    DEBUG(
        std::cout << "========" << std::endl;
    )
    
    //get the right shadows
    for(std::list<Seed>::iterator pSeed = pSeeds->begin(); pSeed != pSeeds->end(); pSeed++)
        vShadows.push_back(getRightShadow(pSeed));

    //perform the line sweep algorithm on the right shadows
    linesweep(vShadows, pSeeds);

    pSeeds->sort(
            [](const Seed& xA, const Seed& xB)
            {
                if(xA.start_ref() == xB.start_ref())
                    return xA.start() < xB.start();
                return xA.start_ref() < xB.start_ref();
            }//lambda
        );//sort function call

    //copy pSeeds since we want to return it
    return pSeeds;
}//function

void exportLinesweep()
{
    //export the LineSweepContainer class
	boost::python::class_<
        LineSweep, 
        boost::python::bases<CppModule>,
        std::shared_ptr<LineSweep>
    >(
        "LineSweep",
        "Uses linesweeping to remove contradicting "
        "matches within one strip of consideration.\n"
        "\n"
        "Execution:\n"
        "	Expects query, ref, strip_vec as input.\n"
        "		query: the query as NucleotideSequence\n"
        "		ref: the reference sequence as Pack\n"
        "		strip_vec: the areas that shall be evaluated as StripOfConsiderationVector\n"
        "	returns strip_vec.\n"
        "		strip_vec: the evaluated areas\n"
    );
	boost::python::implicitly_convertible< 
		std::shared_ptr<LineSweep>,
		std::shared_ptr<CppModule> 
	>();
}//function