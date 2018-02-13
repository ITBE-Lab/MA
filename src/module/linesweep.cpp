#include "module/linesweep.h"
using namespace libMA;

extern int iGap;
extern int iExtend;
extern int iMatch;

ContainerVector LinearLineSweep::getInputType() const
{
    return ContainerVector{
            std::shared_ptr<Container>(new Seeds())
        };
}//function

std::shared_ptr<Container> LinearLineSweep::getOutputType() const
{
    return std::shared_ptr<Container>(new Seeds());
}//function

/**
 * determine the start and end positions this match casts on the left border of the given bucket
 * pxMatch is the container this match is stored in.
 */
LinearLineSweep::ShadowInterval LinearLineSweep::getLeftShadow(Seeds::iterator pSeed) const
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
LinearLineSweep::ShadowInterval LinearLineSweep::getRightShadow(Seeds::iterator pSeed) const
{
    return ShadowInterval(
            pSeed->start_ref(),
            pSeed->end() - (int64_t)pSeed->start_ref(),
            pSeed
        );
}//function




void LinearLineSweep::linesweep(
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
    std::list<ShadowInterval> xItervalEnds = std::list<ShadowInterval>();

    //this is the line sweeping part
    for(ShadowInterval& rInterval : vShadows)
    {
        DEBUG(
            std::cout << "Current Sweep position: " << rInterval.start() << std::endl;
            std::cout << "\tat start of interval " << rInterval.start() <<
                ", " << rInterval.end() << std::endl;
            std::cout << "\t(start_ref, end_ref; start_query, end_query) " 
                << rInterval.pSeed->start_ref() << ", " << rInterval.pSeed->end_ref() << "; "
                << rInterval.pSeed->start() << ", " << rInterval.pSeed->end() << std::endl;
        )
        //check weather there are contradictions to the current seed
        if(!xItervalEnds.empty() && xItervalEnds.front().end() >= rInterval.end())
        {
            //yes there are some!

            // special case: we have a duplicate seed => we do not want to remove both instances
            // cause mereley a duplicate is not a contradiction
            // note: we do not need to check for duplicates later due to the ordering of the seeds
            if(rInterval.within(xItervalEnds.front()))
            {
                xItervalEnds.front().remove(pSeeds);
                continue;
            }//if
            DEBUG(
                std::cout << "\tremoving: " << xItervalEnds.front().end() << std::endl;
            )
            // delete the first seed
            xItervalEnds.front().remove(pSeeds);
            // check if following seeds need also be deleted
            auto iterator = xItervalEnds.begin();
            // we dealed with the first element already
            ++iterator;
            while(iterator != xItervalEnds.end() && iterator->end() >= rInterval.end())
            {
                DEBUG(
                    std::cout << "\tremoving: " << iterator->end() << std::endl;
                )
                iterator->remove(pSeeds);
                iterator = xItervalEnds.erase(iterator);
            }//while
            // remove the current seed since there were contradictions to it.
            rInterval.remove(pSeeds);
        }//if

        //there are no contradictions -> remember the end of the seed
        else
            xItervalEnds.push_front(rInterval);
    }//for
}//function

std::shared_ptr<Container> LinearLineSweep::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
    // copy the input
    // @todo unecessary copy!!! 
    // (well it's necessary since python might delete the old datastructure...)
    std::shared_ptr<Seeds> pSeeds = std::shared_ptr<Seeds>(new Seeds(
        std::static_pointer_cast<Seeds>((*vpInput)[0])));


    std::vector<ShadowInterval> vShadows = {};

    // get the left shadows
    for(Seeds::iterator pSeed = pSeeds->begin(); pSeed != pSeeds->end(); pSeed++)
        vShadows.push_back(getLeftShadow(pSeed));

    // perform the line sweep algorithm on the left shadows
    linesweep(vShadows, pSeeds);

    vShadows.clear();

    DEBUG(
        std::cout << "========" << std::endl;
    )

    // get the right shadows
    for(Seeds::iterator pSeed = pSeeds->begin(); pSeed != pSeeds->end(); pSeed++)
        vShadows.push_back(getRightShadow(pSeed));

    // perform the line sweep algorithm on the right shadows
    linesweep(vShadows, pSeeds);

    pSeeds->bConsistent = true;

    if(pSeeds->size() <= 1)
        return pSeeds;

    // seeds need to be sorted for the following steps
    std::sort(
            pSeeds->begin(), pSeeds->end(),
            [](const Seed& xA, const Seed& xB)
            {
                if(xA.start_ref() == xB.start_ref())
                    return xA.start() < xB.start();
                return xA.start_ref() < xB.start_ref();
            }//lambda
        );//sort function call

    /*
     * FILTER:
     * do a linear gap cost estimation 
     * this does not improve quality but performance
     * since we might remove some seeds that are too far from another
     * 
     * this is much like the backtracking in SW/NW
     * we make sure that we can only have a positive score
     */
    //query position of last encountered seed
    nucSeqIndex uiLastQ = pSeeds->front().start();
    //reference position of last encountered seed
    nucSeqIndex uiLastR = pSeeds->front().start_ref();
    //running score
    nucSeqIndex uiScore = 0;
    //maximal score
    nucSeqIndex uiMaxScore = 0;
    //points to the last seed that had a score of 0
    Seeds::iterator pLastStart = pSeeds->begin();
    //outcome
    Seeds::iterator pOptimalStart = pSeeds->begin();
    //outcome
    Seeds::iterator pOptimalEnd = pSeeds->end();

    /*
     * the goal of this loop is to set pOptimalStart & end correctly
     */
    for(Seeds::iterator pSeed = pSeeds->begin(); pSeed != pSeeds->end(); pSeed++)
    {
        assert(pSeed->start() <= pSeed->end());
        assert(pSeed->size() <= 10000);
        //adjust the score correctly
        uiScore += iMatch * pSeed->getValue();
        /*
         * we need to extract the gap between the seeds in order to do that.
         * we will assume that any gap can be spanned by a single indel
         * while all nucleotides give matches.
         * Therefore we need to get the width x and height y of the rectangular gap
         * number of matches equals min(x, y)
         * gap size equals |x-y|
         */
        nucSeqIndex uiGap = 0;
        if(pSeed->start() > uiLastQ)
            uiGap = pSeed->start() - uiLastQ;
        if(pSeed->start_ref() > uiLastR)
        {
            if(pSeed->start_ref() - uiLastR < uiGap)
            {
                uiGap -= pSeed->start_ref() - uiLastR;
                if(optimisticGapEstimation)
                    uiScore += iMatch * (pSeed->start_ref() - uiLastR);
            }//if
            else
            {
                if(optimisticGapEstimation)
                    uiScore += iMatch * uiGap;
                uiGap = (pSeed->start_ref() - uiLastR) - uiGap;
            }//else
        }//if
        uiGap *= iExtend;
        if(uiGap > 0)
            uiGap += iGap;
        if(uiScore < uiGap)
        {
            uiScore = 0;
            pLastStart = pSeed;
        }//if
        else
            uiScore -= uiGap;
        if(uiScore > uiMaxScore)
        {
            uiMaxScore = uiScore;
            pOptimalStart = pLastStart;
            pOptimalEnd = pSeed;
        }//if
    }//for

    /*
     * we need to increment both iterator but check if they extend past the end of pSeeds
     * (check twice because incrementing an iterator pointing to end will
     * result in undefined behaviour)
     *
     * We then simply remove all known suboptimal seeds
     * this is an optimistic estimation so some suboptimal regions might remain
     * 
     * NOTE: order is significant:
     * """
     * --- Iterator validity ---
     * Iterators, pointers and references pointing to position (or first) and beyond are
     * invalidated, with all iterators, pointers and references to elements before position
     * (or first) are guaranteed to keep referring to the same elements they were referring to
     * before the call.
     * """
     */
    if(pOptimalEnd != pSeeds->end())
        if(++pOptimalEnd != pSeeds->end())
            pSeeds->erase(pOptimalEnd, pSeeds->end());
    if(pOptimalStart != pSeeds->end())
        if(++pOptimalStart != pSeeds->end())
            pSeeds->erase(pSeeds->begin(), pOptimalStart);

    return pSeeds;
}//function

void exportLinesweep()
{
    //export the LineSweepContainer class
    boost::python::class_<
            LinearLineSweep, 
            boost::python::bases<Module>,
            std::shared_ptr<LinearLineSweep>
        >("LinearLineSweep")
        .def_readwrite("optimistic_gap_estimation", &LinearLineSweep::optimisticGapEstimation)
    ;
    boost::python::implicitly_convertible< 
        std::shared_ptr<LinearLineSweep>,
        std::shared_ptr<Module> 
    >();
}//function