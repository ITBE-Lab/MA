#ifndef LINESWEEP
#define LINESWEEP

#include "module.h"
#include "interval.h"
#include "seed.h"
#include "graphicalMethod.h"
#include <memory>
#include "balancedSearchTree.h"
#include "meta_programming.h"

/** each perfect match "casts a shadow" at the left and right border of the strip
 * each shadow is stored in one of these data structures.
*/
class ShadowInterval: public Interval<nucSeqIndex>{
private:
	////while swiping the interfering shadows will get stored in this list
    std::shared_ptr<std::list<ShadowInterval*>> pInterferingIntervals;
    ////the interval this one interfers with
    ShadowInterval* pIInterferWith;
    ////the seed this interval corresponds to (used to delete the seed in case this is necessary)
    std::list<Seed>::iterator pSeed;
	////the total score of the interfering shadows
    unsigned int iScoreInterfering;
public:
    ShadowInterval(
            nucSeqIndex iBegin, 
            nucSeqIndex iSize, 
            std::list<Seed>::iterator pSeed
        )
            :
        Interval(iBegin, iSize),
        pInterferingIntervals(new std::list<ShadowInterval*>()),
        pIInterferWith(),
        pSeed(pSeed),
        iScoreInterfering(0)
    {}//constructor

    ShadowInterval( const ShadowInterval& rOther )
            :
        Interval(rOther),
        pInterferingIntervals(rOther.pInterferingIntervals),
        pIInterferWith(rOther.pIInterferWith),
        pSeed(rOther.pSeed),
        iScoreInterfering(rOther.iScoreInterfering)
    {}//constructor

    nucSeqIndex non_overlap(ShadowInterval rInterval)
    {
        return std::max(
                (*pInterferingIntervals->back())->end_ref() > rInterval->start_ref() ?
                (*pInterferingIntervals->back())->end_ref() - rInterval->start_ref() :
                0
                ,
                (*pInterferingIntervals->back())->end() > rInterval->start() ?
                (*pInterferingIntervals->back())->end() - rInterval->start() :
                0
            );
    }//function
    
    void addInterferingInterval(ShadowInterval& rInterval)
    {
        DEBUG(
            std::cout << "\tis Interfering with interval: " << start() << ", " << end() << std::endl;
        )
        if(pInterferingIntervals->empty())
            iScoreInterfering += rInterval.pSeed->getValue();
        else
            iScoreInterfering += non_overlap(rInterval);
        pInterferingIntervals->push_back(&rInterval);
        
    }//function

    void removeInterferingIntervals(
            std::list<Seed>& rSeeds,
            std::shared_ptr<StripOfConsideration> pStrip
        )
    {
        for(ShadowInterval* pInterval : *pInterferingIntervals)
            pInterval->removeInterferingIntervals(rSeeds, pStrip);
        //reached an already invalidated itterator
        if(pSeed == rSeeds.end())
        {
            DEBUG(
                std::cout << "\treached invalid itterator" << std::endl;
            )
            return;
        }
        DEBUG(
            std::cout << "\terasing: " << start() << ", " << end();
            std::cout << " seed: " << pSeed->start() << ", " << pSeed->end() << std::endl;
        )
        pStrip->subtractFromValue(pSeed->getValue());
        rSeeds.erase(pSeed);
        pSeed = rSeeds.end();
        pInterferingIntervals->clear();
    }//function

    void removeSeedIfNecessary(
            std::list<Seed>& rSeeds, 
            std::shared_ptr<StripOfConsideration> pStrip
        )
    {
        //reached an already invalidated itterator
        if(pSeed == rSeeds.end())
        {
            DEBUG(
                std::cout << "\treached invalid itterator" << std::endl;
            )
            return;
        }
        //the many interfering intervals are more valuable => remove this interval
        if(pSeed->getValue() < iScoreInterfering)
        {
            //update the score of the interval outside this one
            if(pIInterferWith != nullptr)
                pIInterferWith->iScoreInterfering += iScoreInterfering - pSeed->getValue();
            pStrip->subtractFromValue(pSeed->getValue());
            //remove this seed
            rSeeds.erase(pSeed);
            pSeed = rSeeds.end();
            DEBUG(
                std::cout << "\tis less valuable than interfering ones" << std::endl;
            )
        }//if
        //this interval is more valuable than the interfering ones => remove the interfering ones
        else
        {
            DEBUG(
                std::cout << "\tis more valuable than interfering ones" << std::endl;
            )
            for(ShadowInterval* pInterval : *pInterferingIntervals)
                pInterval->removeInterferingIntervals(rSeeds, pStrip);
            pInterferingIntervals->clear();
        }//else
    }//function

    inline const std::list<Seed>::iterator operator*() const
    {
        return pSeed;
    }//operator
    inline const std::list<Seed>::iterator operator->() const
    {
        return pSeed;
    }//operator

    inline ShadowInterval& operator=(const ShadowInterval& rxOther)
    {
        Interval::operator=(rxOther);
        pInterferingIntervals = rxOther.pInterferingIntervals;
        pIInterferWith = rxOther.pIInterferWith;
        pSeed = rxOther.pSeed;
        return *this;
    }// operator

    //required for the search tree
    inline bool operator<=(const ShadowInterval rOther) const
    {
        return end() <= rOther.end();
    }//operator

    //required for the search tree
    inline bool operator<(const ShadowInterval rOther) const
    {
        return end() < rOther.end();
    }//operator

    //required for the search tree
    inline bool operator>(const ShadowInterval rOther) const
    {
        return end() > rOther.end();
    }//operator
};//class

class LineSweep: public Module
{
private:
    void linesweep(
            std::vector<ShadowInterval>& vShadows, 
            std::list<Seed>& rSeeds,
            std::shared_ptr<StripOfConsideration> pStrip
        );

    ShadowInterval getRightShadow(
            nucSeqIndex iBucketStart,
            std::list<Seed>::iterator pSeed,
            nucSeqIndex iQueryLength
        ) const;

    ShadowInterval getLeftShadow(
            nucSeqIndex uiBucketStart,
            std::list<Seed>::iterator pSeed,
            nucSeqIndex iRefSize
        ) const;
public:

	LineSweep(){}//constructor

	std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> pInput);

    std::vector<ContainerType> getInputType();

    std::vector<ContainerType> getOutputType();
};//class

void exportLinesweep();

#endif