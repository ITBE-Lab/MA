#ifndef LINESWEEP
#define LINESWEEP

#include "module.h"
#include "interval.h"
#include "seed.h"
#include "graphicalMethod.h"
#include <memory>
#include "balancedSearchTree.h"

/** each perfect match "casts a shadow" at the left and right border of the strip
 * each shadow is stored in one of these data structures.
*/
class ShadowInterval: public Interval<nucSeqIndex>{
private:
	////while swiping the interfering shadows will get stored in this list
    std::list<ShadowInterval*> lInterferingIntervals;
    ////the interval this one interfers with
    ShadowInterval* pIInterferWith;
    ////the seed this interval corresponds to (used to delete the seed in case this is necessary)
    std::list<Seed>::iterator pSeed;
	////the total score of the interfering shadows
    unsigned int iScoreInterfering;
public:
    ShadowInterval(
            nucSeqIndex iBegin, 
            nucSeqIndex iEnd, 
            std::list<Seed>::iterator pSeed
        )
            :
        Interval(iBegin, iEnd),
        lInterferingIntervals(),
        pIInterferWith(),
        pSeed(pSeed),
        iScoreInterfering(0)
    {}//constructor

    nucSeqIndex non_overlap(ShadowInterval* pInterval)
    {
        return std::max(
                lInterferingIntervals.back()->pSeed->end_ref() > pInterval->pSeed->start_ref() ?
                lInterferingIntervals.back()->pSeed->end_ref() - pInterval->pSeed->start_ref() :
                0
                ,
                lInterferingIntervals.back()->pSeed->end() > pInterval->pSeed->start() ?
                lInterferingIntervals.back()->pSeed->end() - pInterval->pSeed->start() :
                0
            );
    }//function
    
    void addInterferingInterval(ShadowInterval* pInterval)
    {
        if(lInterferingIntervals.empty())
            iScoreInterfering += pInterval->pSeed->getValue();
        else
            iScoreInterfering += non_overlap(pInterval);
        lInterferingIntervals.push_back(pInterval);
    }//function

    void removeInterferingIntervals(std::list<Seed>& rSeeds)
    {
        for(ShadowInterval* pInterval : lInterferingIntervals)
            pInterval->removeInterferingIntervals(rSeeds);
        rSeeds.erase(pSeed);
    }//function

    void removeSeedIfNecessary(std::list<Seed>& rSeeds)
    {
        //the many interfering intervals are more valuable => remove this interval
        if(pSeed->getValue() < iScoreInterfering)
        {
            //update the score of the interval outside this one
            if(pIInterferWith != nullptr)
                pIInterferWith->iScoreInterfering += iScoreInterfering - pSeed->getValue();
            //remove this seed
            rSeeds.erase(pSeed);
        }//if
        //this interval is more valuable than the interfering ones => remove the interfering ones
        else
            removeInterferingIntervals(rSeeds);
    }//function

    inline const Seed& operator*() const
    {
        return *pSeed;
    }//operator
    inline const Seed& operator->() const
    {
        return *pSeed;
    }//operator

    inline ShadowInterval& operator=(const ShadowInterval& rxOther)
    {
        Interval::operator=(rxOther);
        lInterferingIntervals = rxOther.lInterferingIntervals;
        pSeed = rxOther.pSeed;
        return *this;
    }// operator

    //required for the search tree
    inline bool operator<(const ShadowInterval* pOther) const
    {
        return end() < pOther->end();
    }//operator

    //required for the search tree
    inline bool operator<(const unsigned int value) const
    {
        return end() < value;
    }//operator
};//class

class LineSweep: public Module
{
private:
    void linesweep(std::vector<ShadowInterval>& vShadows, std::list<Seed>& rSeeds);

    ShadowInterval getRightShadow(
            nucSeqIndex iBucketStart,
            std::list<Seed>::iterator pSeed,
            nucSeqIndex iBucketSize,
            nucSeqIndex iQueryLength
        ) const;

    ShadowInterval getLeftShadow(
            nucSeqIndex uiBucketStart,
            std::list<Seed>::iterator pSeed,
            nucSeqIndex uiBucketSize,
            nucSeqIndex uiQueryLength
        ) const;
public:

	LineSweep(){}//constructor

	std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> pInput);

    std::vector<ContainerType> getInputType();

    std::vector<ContainerType> getOutputType();
};//class

void exportLinesweep();

#endif