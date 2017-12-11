/** 
 * @file linesweep.h
 * @brief Implements the linesweep @ref Module "module"
 * @author Markus Schmidt
 */
#ifndef LINESWEEP
#define LINESWEEP

#include "segmentList.h"
#include "searchTree.h"
#include "module.h"

namespace libLAuS
{
    // deprecated old linesweep code
    #if 0
    /**
     * @brief The shadow of a Seed.
     * @details
     * Each perfect match "casts a shadow" at the left and right border of the strip.
     * Each shadow is stored in one of these data structures.
     * @note DEPRECATED!
     */
    class ShadowInterval: public Interval<int64_t>{
    public:
        /// @brief While swiping the interfering shadows will get stored in this list.
        std::shared_ptr<std::list<ShadowInterval*>> pInterferingIntervals;
        /// @brief old complicated algorithm required that
        std::shared_ptr<std::list<ShadowInterval*>> pInterferingIntervals2ndOrder;
        /// @brief The interval this one interferes with.
        ShadowInterval* pIInterferWith;
        /// @brief old complicated algorithm required that
        std::list<std::tuple<ShadowInterval*, unsigned int>> lIInterferWith2ndOrder;
        /// @brief The seed this interval corresponds to (used to delete the seed in case this is necessary).
        std::list<Seed>::iterator pSeed;
        /// @brief The total score of the interfering shadows.
        unsigned int iScoreInterfering;
        /// @brief old complicated algorithm required that
        unsigned int iScoreInterfering2ndOrder;
        /// @brief The score added by this interval to the outer interfering one
        unsigned int iInterferingSelf;

        /**
         * @brief Creates a new shadow.
         * @details
         * The linesweep algorithm disables seeds.
         * Therefore the iterator is required in order to delete the respective seed from its list.
         */
        ShadowInterval(
                int64_t iBegin, 
                int64_t iSize, 
                std::list<Seed>::iterator pSeed
            )
                :
            Interval(iBegin, iSize),
            pInterferingIntervals(new std::list<ShadowInterval*>()),
            pInterferingIntervals2ndOrder(new std::list<ShadowInterval*>()),
            pIInterferWith(),
            lIInterferWith2ndOrder(),
            pSeed(pSeed),
            iScoreInterfering(0),
            iScoreInterfering2ndOrder(0),
            iInterferingSelf(0)
        {}//constructor

        /**
         * @brief Copy constructor
         */
        ShadowInterval( const ShadowInterval& rOther )
                :
            Interval(rOther),
            pInterferingIntervals(rOther.pInterferingIntervals),
            pInterferingIntervals2ndOrder(rOther.pInterferingIntervals2ndOrder),
            pIInterferWith(rOther.pIInterferWith),
            lIInterferWith2ndOrder(rOther.lIInterferWith2ndOrder),
            pSeed(rOther.pSeed),
            iScoreInterfering(rOther.iScoreInterfering),
            iScoreInterfering2ndOrder(rOther.iScoreInterfering2ndOrder),
            iInterferingSelf(rOther.iInterferingSelf)
        {}//copy constructor

        /**
         * @brief Returns the length of this interval that is not overlapped by the given one.
         */
        int64_t non_overlap(ShadowInterval rInterval)
        {
            if(pInterferingIntervals->empty())
                return size();
            return std::max(
                    pInterferingIntervals->back()->pSeed->end_ref() > rInterval.pSeed->start_ref() ?
                    pInterferingIntervals->back()->pSeed->end_ref() - rInterval.pSeed->start_ref() :
                    0
                    ,
                    pInterferingIntervals->back()->pSeed->end() > rInterval.pSeed->start() ?
                    pInterferingIntervals->back()->pSeed->end() - rInterval.pSeed->start() :
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
                rInterval.iInterferingSelf = rInterval.pSeed->getValue();
            else
                rInterval.iInterferingSelf = non_overlap(rInterval);
            iScoreInterfering += rInterval.iInterferingSelf;
            pInterferingIntervals->push_back(&rInterval);
            rInterval.pIInterferWith = this;
        }//function

        
        void add2ndOrderInterferingInterval(ShadowInterval& rInterval)
        {
            DEBUG(
                std::cout << "\tis Interfering with interval (2nd order): " << start() << ", ";
                std::cout << end() << std::endl;
            )
            unsigned int adding = 0;
            if(pInterferingIntervals2ndOrder->empty())
                adding = rInterval.pSeed->getValue();
            else
                adding = non_overlap(rInterval);
            iScoreInterfering2ndOrder += adding;
            pInterferingIntervals2ndOrder->push_back(&rInterval);
            rInterval.lIInterferWith2ndOrder.push_back(std::make_tuple(this, adding));
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
            assert(pInterferingIntervals != nullptr);
            for(ShadowInterval* pInterval : *pInterferingIntervals)
            {
                assert(pInterval != nullptr);
                assert(pInterval != this);
                pInterval->removeInterferingIntervals(pSeeds);
            }//for
            //reached an already invalidated iterator
            if(pSeed == pSeeds->end())
            {
                DEBUG(
                    std::cout << "\treached invalid iterator" << std::endl;
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
            for(std::tuple<ShadowInterval*, unsigned int> tup : lIInterferWith2ndOrder)
            {
                std::get<0>(tup)->iScoreInterfering2ndOrder -= std::get<1>(tup);
            }//for
            lIInterferWith2ndOrder.clear();
        }//function

        
        void removeInterferingIntervals2ndOrder(
                std::shared_ptr<Seeds> pSeeds
            )
        {
            for(ShadowInterval* pInterval : *pInterferingIntervals2ndOrder)
                pInterval->removeInterferingIntervals2ndOrder(pSeeds);
            //reached an already invalidated iterator
            if(pSeed == pSeeds->end())
            {
                DEBUG(
                    std::cout << "\treached invalid iterator (2nd order)" << std::endl;
                )
                return;
            }
            DEBUG(
                std::cout << "\terasing (2nd order): " << start() << ", " << end();
                std::cout << " seed (2nd order): " << pSeed->start() << ", " << pSeed->end() << std::endl;
            )
            pSeeds->erase(pSeed);
            pSeed = pSeeds->end();
            pInterferingIntervals2ndOrder->clear();
            //adjust scores for outer intervals
            for(std::tuple<ShadowInterval*, unsigned int> tup : lIInterferWith2ndOrder)
            {
                std::get<0>(tup)->iScoreInterfering2ndOrder -= std::get<1>(tup);
            }//for
            //adjust 2nd order score
            if(pIInterferWith != nullptr && pIInterferWith->pIInterferWith != nullptr)
                //pIInterferWith MUST be deactivated already!
                //but if it interfers with another interval we have to update the score of the outer one
                pIInterferWith->pIInterferWith->iScoreInterfering -= iInterferingSelf;
            lIInterferWith2ndOrder.clear();
        }//function


        /**
         * @brief Removes this or all interfering seeds.
         * @details
         * Checks weather this seed is more valuable than the interfering ones.
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
            if(pSeed->getValue() < iScoreInterfering + iScoreInterfering2ndOrder)
            {
                //update the score of the interval outside this one
                if(pIInterferWith != nullptr)
                    pIInterferWith->iScoreInterfering += iScoreInterfering - iInterferingSelf;
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
                    std::cout << "\tremoving " << pInterferingIntervals->size() << " + ";
                    std::cout << pInterferingIntervals2ndOrder->size() << " intervals." << std::endl;
                )
                for(ShadowInterval* pInterval : *pInterferingIntervals)
                    pInterval->removeInterferingIntervals(pSeeds);
                pInterferingIntervals->clear();
                for(ShadowInterval* pInterval : *pInterferingIntervals2ndOrder)
                    pInterval->removeInterferingIntervals2ndOrder(pSeeds);
                pInterferingIntervals2ndOrder->clear();
            }//else
        }//function

        inline const ShadowInterval* getIInterferWith() const
        {
            return pIInterferWith;
        }//function


        
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
        inline bool operator<=(const ShadowInterval& rOther) const
        {
            return end() <= rOther.end();
        }//operator

        /**
         * @brief Compares the Interval ends.
         * @details
         * Required for the search tree.
         */
        inline bool operator<(const ShadowInterval& rOther) const
        {
            return end() < rOther.end();
        }//operator

        /**
         * @brief Compares the Interval ends.
         * @details
         * Required for the search tree.
         */
        inline bool operator>(const ShadowInterval& rOther) const
        {
            return end() > rOther.end();
        }//operator
    };//class

    /**
     * @details
     * small class in order to have a operator< on a pointer
     */
    class ShadowIntervalPtr
    {
    public:
        ShadowInterval *p;

        ShadowIntervalPtr(ShadowInterval* p)
                :
            p(p)
        {}//constructor
        
        inline bool operator<(const ShadowIntervalPtr& rOther) const
        {
            return *p < *rOther.p;
        }//operator
        
        inline bool operator==(const ShadowIntervalPtr& rOther) const
        {
            return *p == *rOther.p;
        }//operator

        /**
         * @brief Provides easy access to the seed within this shadow.
         */
        inline ShadowInterval* operator->() const
        {
            return p;
        }//operator
    };//class

    /**
     * @brief Implements the linesweep algorithm.
     * @ingroup module
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
        ShadowInterval getLeftShadow(std::list<Seed>::iterator pSeed) const;

        /**
        * @brief Returns the left shadow of a seed.
        * @details
        * "Casts" the shadow against the border of the considered area.
        */
        ShadowInterval getRightShadow(std::list<Seed>::iterator pSeed) const;
    public:

        //overload
        std::shared_ptr<Container> execute(ContainerVector pInput);

        //overload
        ContainerVector getInputType() const;

        //overload
        std::shared_ptr<Container> getOutputType() const;

        std::string getName() const
        {
            return "LineSweep";
        }
    };//class
    #endif


    /**
     * @brief The shadow of a Seed.
     * @details
     * Each perfect match "casts a shadow" at the left and right border of the strip.
     * Each shadow is stored in one of these data structures.
     */
    class ShadowInterval2: public Interval<int64_t>{
    public:
        Seeds::iterator pSeed;

        /**
         * @brief Creates a new shadow.
         * @details
         * The linesweep algorithm disables seeds.
         * Therefore the iterator is required in order to delete the respective seed from its list.
         */
        ShadowInterval2(
                int64_t iBegin, 
                int64_t iSize, 
                Seeds::iterator pSeed
            )
                :
            Interval(iBegin, iSize),
            pSeed(pSeed)
        {}//constructor

        /**
         * @brief Copy constructor
         */
        ShadowInterval2( const ShadowInterval2& rOther )
                :
            Interval(rOther),
            pSeed(rOther.pSeed)
        {}//copy constructor

        /**
         * @brief Removes this or all interfering seeds.
         * @details
         * Checks weather this seed is more valuable than the interfering ones.
         */
        void remove(
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
            pSeeds->erase(pSeed);
            pSeed = pSeeds->end();
        }//function

        bool within(const ShadowInterval2& rOther)
        {
            return start() >= rOther.start() && end() <= rOther.end();
        }//function
    };//class


    /**
     * @brief Implements the LinearLineSweep algorithm.
     * @ingroup module
     * @details
     * Removes all contradicting seeds.
     * This should only be used in combination with the StripOfConsideration module.
     */
    class LinearLineSweep: public Module
    {
    private:
        /**
        * @brief Implements the linesweep algorithm.
        * @details
        * The algorithm has to be run on left and right shadows,
        * therefore it is provided as individual function.
        */
        void linesweep(
                std::vector<ShadowInterval2>& vShadows, 
                std::shared_ptr<Seeds> pSeeds
            );

        /**
        * @brief Returns the left shadow of a seed.
        * @details
        * "Casts" the left shadows.
        */
        ShadowInterval2 getLeftShadow(Seeds::iterator pSeed) const;

        /**
        * @brief Returns the right shadow of a seed.
        * @details
        * "Casts" the right shadows.
        */
        ShadowInterval2 getRightShadow(Seeds::iterator pSeed) const;
    public:

        //overload
        std::shared_ptr<Container> execute(ContainerVector pInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - Seeds
         */
        ContainerVector getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - Seeds
         */
        std::shared_ptr<Container> getOutputType() const;

        std::string getName() const
        {
            return "LineSweep2";
        }
    };//class
}//namespace libLAuS

/**
 * @brief Exposes the LineSweep @ref Module "module" to boost python.
 * @ingroup export
 */
void exportLinesweep();

#endif