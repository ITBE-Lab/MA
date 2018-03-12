/** 
 * @file segment.h
 * @brief Implements a the IntervalTree used for segmentation and various other related classes.
 * @author Markus Schmidt
 */
#ifndef INTERVALTREE_H
#define INTERVALTREE_H

#include "container/fMIndex.h"
#include <thread>
#include "container/seed.h"

#define confMETA_MEASURE_DURATION ( 1 )

namespace libMA
{

    /**
     * @brief A Suffix Array Segment.
     * @details
     * A Suffix Array Segment is made up of two Intervals.
     * @li @c a SAInterval.
     * @li @c a Interval representing the position of the sequence on the query.
     * @ingroup container
     */
    class Segment: public Container, public Interval<nucSeqIndex> {
    public:
        SAInterval xSaInterval;
        /**
        * @brief Creates a new Segment.
        * @details Creates a new Segment on the base of a SAInterval and the 
        * respective indices on the quey.
        */
        Segment(nucSeqIndex uiStart, nucSeqIndex uiSize, SAInterval xSaInterval)
                :
            Interval(uiStart, uiSize),
            xSaInterval(xSaInterval)
        {}//constructor

        Segment()
                :
            Interval(),
            xSaInterval()
        {}//constructor

        Segment(const Segment& other)
                :
            Interval(other),
            xSaInterval(other.xSaInterval)
        {}//copy constructor


        //overload
        bool canCast(std::shared_ptr<Container> c) const
        {
            return std::dynamic_pointer_cast<Segment>(c) != nullptr;
        }//function

        //overload
        std::string getTypeName() const
        {
            return "Segment";
        }//function

        //overload
        std::shared_ptr<Container> getType() const
        {
            return std::shared_ptr<Container>(new Segment());
        }//function

        /**
         * @brief The bwt interval within.
         * @returns the bwt interval within.
         */
        const SAInterval& saInterval() const
        {
            return xSaInterval;
        }//function
        
        /**
         * @brief Copys from another Segment.
         * @note override
         */
        inline Segment& operator=(const Segment& rxOther)
        {
            Interval::operator=(rxOther);
            xSaInterval = rxOther.xSaInterval;
            return *this;
        }// operator

    }; // class ( Segment )


    /**
     * @brief The segment tree.
     * @details
     * The segment "tree" is actually a doubly linked list.
     * The tree only exists logically,
     * meaning that the segments within the list represent the first layer of the tree initally.
     * Then after each iteration, the segments within the list represent the next layer down of the 
     * tree.
     * @ingroup container
     */
    class SegmentVector : public std::vector<std::shared_ptr<Segment>>, public Container{

    public:

        template< class InputIt >
        SegmentVector(InputIt xBegin, InputIt xEnd
            )
            :
            vector(xBegin, xEnd)
        {}//iterator constructor

        SegmentVector()
        {}//default constructor

        //overload
        bool canCast(std::shared_ptr<Container> c) const
        {
            return std::dynamic_pointer_cast<SegmentVector>(c) != nullptr;
        }//function

        //overload
        std::string getTypeName() const
        {
            return "SegmentVector";
        }//function

        //overload
        std::shared_ptr<Container> getType() const
        {
            return std::shared_ptr<SegmentVector>(new SegmentVector());
        }//function

        /**
         * @brief Extracts all seeds from the tree.
         * @details
         * Calls fDo for all recorded hits.
         * fDo shall return false to terminate the iteration shall end.
         * @Note pushBackBwtInterval records an interval of hits
         */
        void forEachSeed(
                std::shared_ptr<FMIndex> pxFMIndex,
                unsigned int uiMAxAmbiguity,
                bool bSkip,
                std::function<bool(Seed s)> fDo
            )
        {
            //iterate over all the intervals that have been recorded using pushBackBwtInterval()
            for (std::shared_ptr<Segment> pSegment : *this)
            {
                if(pSegment == nullptr)
                {
                    DEBUG(
                        std::cout << "WARNING: found nullptr in SegmentVector!" << std::endl;
                    )
                    continue;
                }
                //if the interval contains more than uiMAxAmbiguity hits it's of no importance and will produce nothing but noise

                //if bSkip is not set uiJump by is used to not return more than uiMAxAmbiguity
                
                t_bwtIndex uiJumpBy = 1;
                if (pSegment->saInterval().size() > uiMAxAmbiguity && uiMAxAmbiguity != 0)
                {
                    if (bSkip)
                        continue;
                    uiJumpBy = pSegment->saInterval().size() / uiMAxAmbiguity; 
                }//if

                //iterate over the interval in the BWT
                for (
                        auto ulCurrPos = pSegment->saInterval().start(); 
                        ulCurrPos < pSegment->saInterval().end(); 
                        ulCurrPos += uiJumpBy
                    )
                {
                    //calculate the referenceIndex using pxUsedFmIndex->bwt_sa() and call fDo for every match individually
                    nucSeqIndex ulIndexOnRefSeq = pxFMIndex->bwt_sa(ulCurrPos);
                    //call the given function
                    if(!fDo(Seed(
                            pSegment->start(),
                            pSegment->size() + 1,
                            ulIndexOnRefSeq,
                            pSegment->saInterval().size()
                        )))
                        return;
                }//for
            }//for
        }//function

        /**
         * @brief Extracts all seeds from the segment list.
         */
        std::shared_ptr<Seeds> extractSeeds(
                std::shared_ptr<FMIndex> pxFMIndex, 
                unsigned int uiMAxAmbiguity,
                bool bSkip = true
            )
        {
            std::shared_ptr<Seeds> pRet = std::shared_ptr<Seeds>(new Seeds());
            forEachSeed(
                pxFMIndex,
                uiMAxAmbiguity,
                bSkip,
                [&pRet]
                (Seed s)
                {
                    pRet->push_back(s);
                    return true;
                }//lambda
            );//for each
            return pRet;
        }//function

        
        /**
         * @brief returns the number of seeds
         */
        unsigned int numSeeds(std::shared_ptr<FMIndex> pxFMIndex, unsigned int max_size)
        {
            unsigned int uiTotal = 0;
            for (std::shared_ptr<Segment> pSegment : *this)
                if(max_size == 0 || pSegment->xSaInterval.size() <= max_size)
                    uiTotal += pSegment->xSaInterval.size();
            return uiTotal;
        }// function
    };// class
}// namespace libMA

/**
 * @brief Exposes the SegmentVector to boost python.
 * @ingroup export
 */
void exportIntervalTree();


#endif