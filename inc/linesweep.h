/** 
 * @file linesweep.h
 * @brief Implements the linesweep module
 * @author Markus Schmidt
 */
#ifndef LINESWEEP
#define LINESWEEP

#include "balancedSearchTree.h"
#include "module.h"

/**
 * @brief The shadow of a Seed.
 * @details
 * Each perfect match "casts a shadow" at the left and right border of the strip.
 * Each shadow is stored in one of these data structures.
 */
class ShadowInterval: public Interval<nucSeqIndex>{
private:
	/// @brief While swiping the interfering shadows will get stored in this list.
    std::shared_ptr<std::list<ShadowInterval*>> pInterferingIntervals;
    /// @brief The interval this one interferes with.
    ShadowInterval* pIInterferWith;
    /// @brief The seed this interval corresponds to (used to delete the seed in case this is necessary).
    std::list<Seed>::iterator pSeed;
	/// @brief The total score of the interfering shadows.
    unsigned int iScoreInterfering;
public:
    /**
     * @brief Creates a new shadow.
     * @details
     * The linesweep algorithm disables seeds.
     * Therefore the iterator is required in order to delete the respective seed from its list.
     */
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

    /**
     * @brief Copy constructor
     */
    ShadowInterval( const ShadowInterval& rOther )
            :
        Interval(rOther),
        pInterferingIntervals(rOther.pInterferingIntervals),
        pIInterferWith(rOther.pIInterferWith),
        pSeed(rOther.pSeed),
        iScoreInterfering(rOther.iScoreInterfering)
    {}//copy constructor

    /**
     * @brief Returns the length of this interval that is not overlapped by the given one.
     */
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
    
    /**
     * @brief Records a Interval that interferes with this interval.
     * @details
     * Since the value of all interfering intervals is stored within this interval the score needs
     * to be updated.
     */
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

    /**
     * @brief Removes all intervals that are registered to be interfering with this interval.
     * @details
     * This process is recursive.
     */
    void removeInterferingIntervals(
            std::shared_ptr<Seeds> pSeeds
        )
    {
        for(ShadowInterval* pInterval : *pInterferingIntervals)
            pInterval->removeInterferingIntervals(pSeeds);
        //reached an already invalidated iterator
        if(pSeed == pSeeds->end())
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
        pSeeds->erase(pSeed);
        pSeed = pSeeds->end();
        pInterferingIntervals->clear();
    }//function

    /**
     * @brief Removes this or all interfering seeds.
     * @details
     * Checks weather this seed is more valuable than the interfering ones.
     * Then removes the less valuable ones.
     */
    void removeSeedIfNecessary(
            std::shared_ptr<Seeds> pSeeds
        )
    {
        //reached an already invalidated iterator
        if(pSeed == pSeeds->end())
        {
            DEBUG(
                std::cout << "\treached invalid iterator" << std::endl;
            )
            return;
        }
        //the many interfering intervals are more valuable => remove this interval
        if(pSeed->getValue() < iScoreInterfering)
        {
            //update the score of the interval outside this one
            if(pIInterferWith != nullptr)
                pIInterferWith->iScoreInterfering += iScoreInterfering - pSeed->getValue();
            //remove this seed
            pSeeds->erase(pSeed);
            pSeed = pSeeds->end();
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
                pInterval->removeInterferingIntervals(pSeeds);
            pInterferingIntervals->clear();
        }//else
    }//function

    /**
     * @brief Provides easy access to the seed within this shadow.
     */
    inline const std::list<Seed>::iterator operator*() const
    {
        return pSeed;
    }//operator

    /**
     * @brief Provides easy access to the seed within this shadow.
     */
    inline const std::list<Seed>::iterator operator->() const
    {
        return pSeed;
    }//operator

    
    /**
     * @brief Copys everything from another shadow.
     */
    inline ShadowInterval& operator=(const ShadowInterval& rxOther)
    {
        Interval::operator=(rxOther);
        pInterferingIntervals = rxOther.pInterferingIntervals;
        pIInterferWith = rxOther.pIInterferWith;
        pSeed = rxOther.pSeed;
        return *this;
    }// operator

    /**
     * @brief Compares the Interval ends.
     * @details
     * Required for the search tree.
     */
    inline bool operator<=(const ShadowInterval rOther) const
    {
        return end() <= rOther.end();
    }//operator

    /**
     * @brief Compares the Interval ends.
     * @details
     * Required for the search tree.
     */
    inline bool operator<(const ShadowInterval rOther) const
    {
        return end() < rOther.end();
    }//operator

    /**
     * @brief Compares the Interval ends.
     * @details
     * Required for the search tree.
     */
    inline bool operator>(const ShadowInterval rOther) const
    {
        return end() > rOther.end();
    }//operator
};//class

/**
 * @brief Implements the linesweep algorithm.
 */
class LineSweep: public Module
{
private:
    /**
    * @brief Implements the linesweep algorithm.
    * @details
    * The algorithm has to be run on left and right shadows,
    * therefore it is provided as individual function.
    */
    void linesweep(
            std::vector<ShadowInterval>& vShadows, 
            std::shared_ptr<Seeds> pSeeds
        );

    /**
    * @brief Returns the right shadow of a seed.
    * @details
    * "Casts" the shadow against the border of the considered area.
    */
    ShadowInterval getLeftShadow(
            std::list<Seed>::iterator pSeed,
            nucSeqIndex uiQueryLength
        ) const;

    /**
    * @brief Returns the left shadow of a seed.
    * @details
    * "Casts" the shadow against the border of the considered area.
    */
    ShadowInterval getRightShadow(
            std::list<Seed>::iterator pSeed,
            nucSeqIndex iRefSize
        ) const;
public:

    //overload
	std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> pInput);

    //overload
    std::vector<ContainerType> getInputType();

    //overload
    std::vector<ContainerType> getOutputType();
};//class

/**
 * @brief Exposes the LineSweep Module to boost python.
 */
void exportLinesweep();

#endif