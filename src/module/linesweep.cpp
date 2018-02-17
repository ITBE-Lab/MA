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

std::shared_ptr<std::vector<std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>>>
    LinearLineSweep::linesweep(
        std::shared_ptr<std::vector<
            std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>
        >> pShadows
    )
{
    //sort shadows (increasingly) by start coordinate of the match
    std::sort(
            pShadows->begin(),
            pShadows->end(),
            [](std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex> xA,
                std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex> xB)
            {
                /*
                * sort by the interval starts
                * if two intervals start at the same point the larger one shall be treated first
                */
                if(std::get<1>(xA) == std::get<1>(xB))
                    return std::get<2>(xA) > std::get<2>(xB);
                return std::get<1>(xA) < std::get<1>(xB);
            }//lambda
        );//sort function call

    auto pItervalEnds = std::make_shared<
            std::vector<std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>>
        >();

    nucSeqIndex x = 0;
    //this is the line sweeping part
    for(std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>& xTup : *pShadows)
    {
        if(x < std::get<2>(xTup))
        {
            pItervalEnds->push_back(xTup);
            x = std::get<2>(xTup);
        }//if
        /*
         * sometimes a seeding algorithm may deliver equal seeds
         * in that case these seeds are obviously not contradicting
         *
         * This is not mentioned in the paper since we talk of seed SETs there.
         * A set cannot have the same element twice,
         * however the vector we use in practice can...
         *
         * Anyways we just need to add one simple check for equality here
         * operator!= is not implemented so we rely on ! of ==.
         */
        else if(
                ! pItervalEnds->empty() 
                && ! (*std::get<0>(xTup) == *std::get<0>(pItervalEnds->back())) 
            )
            while(!pItervalEnds->empty() && std::get<2>(pItervalEnds->back()) >= std::get<2>(xTup))
                pItervalEnds->pop_back();
    }//for
    return pItervalEnds;
}//function

std::shared_ptr<Container> LinearLineSweep::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
    std::shared_ptr<Seeds> pSeedsIn = std::static_pointer_cast<Seeds>((*vpInput)[0]);

    auto pShadows = std::make_shared<
            std::vector<std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>>
        >();

    // get the left shadows
    for(Seeds::iterator pSeed = pSeedsIn->begin(); pSeed != pSeedsIn->end(); pSeed++)
        pShadows->push_back(std::make_tuple(
                pSeed,
                pSeed->start(),
                pSeed->end_ref()
            ));

    // perform the line sweep algorithm on the left shadows
    auto pShadows2 = linesweep(pShadows);
    pShadows->clear();

    // get the right shadows
    for(auto &xT : *pShadows2)
        pShadows->push_back(std::make_tuple(
                std::get<0>(xT),
                std::get<0>(xT)->start_ref(),
                std::get<0>(xT)->end()
            ));

    // perform the line sweep algorithm on the right shadows
    pShadows = linesweep(pShadows);

    auto pSeeds = std::make_shared<Seeds>();
    for(auto &xT : *pShadows)
        pSeeds->push_back(*std::get<0>(xT));

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