#ifndef SOC_H
#define SOC_H

#include "container/seed.h"
#include <algorithm>

namespace libMA
{
    /**
     * @brief used to determine more complex orders of SoCs 
     */
    class SoCOrder
    {
    public:
        nucSeqIndex uiAccumulativeLength = 0;
        unsigned int uiSeedAmbiguity = 0;
        unsigned int uiSeedAmount = 0;

        inline void operator+=(const Seed& rS)
        {
            uiSeedAmbiguity += rS.uiAmbiguity;
            uiSeedAmount++;
            uiAccumulativeLength += rS.getValue();
        }//operator

        inline void operator-=(const Seed& rS)
        {
            assert(uiSeedAmbiguity >= rS.uiAmbiguity);
            uiSeedAmbiguity -= rS.uiAmbiguity;
            assert(uiAccumulativeLength >= rS.getValue());
            uiAccumulativeLength -= rS.getValue();
            uiSeedAmount--;
        }//operator

        inline bool operator<(const SoCOrder& rOther) const
        {
            if(uiAccumulativeLength == rOther.uiAccumulativeLength)
                return uiSeedAmbiguity > rOther.uiSeedAmbiguity;
            return uiAccumulativeLength < rOther.uiAccumulativeLength;
        }//operator

        inline void operator=(const SoCOrder& rOther)
        {
            uiAccumulativeLength = rOther.uiAccumulativeLength;
            uiSeedAmbiguity = rOther.uiSeedAmbiguity;
        }//operator
    }; //class


    class SoCPriorityQueue: public Container
    {
    public:
        DEBUG(
            // confirms the respective functions are always called in the correct mode
            bInPriorityMode = false;
        )// DEBUG

        // positions to remember the maxima
        std::vector<std::pair<SoCOrder, std::vector<Seed>::iterator>> vMaxima;
        // last SoC end
        nucSeqIndex uiLastEnd = 0;

        //overload
        inline bool canCast(const std::shared_ptr<Container>& c) const
        {
            return std::dynamic_pointer_cast<SoCPriorityQueue>(c) != nullptr;
        }//function

        //overload
        inline std::string getTypeName() const
        {
            return "SoCPriorityQueue";
        }//function

        //overload
        inline std::shared_ptr<Container> getType() const
        {
            return std::shared_ptr<Container>(new SoCPriorityQueue());
        }//function

        inline bool empty()
        {
            return vMaxima.empty();
        }// function

        inline std::shared_ptr<Seeds> pop()
        {
            DEBUG(assert(bInPriorityMode);)
            //the strip that shall be collected
            std::shared_ptr<Seeds> pSeeds(new Seeds());
            // get the expected amount of seeds in the SoC from the order class and reserve memory
            pSeeds->reserve(vMaxima.front().first.uiSeedAmount);
            //iterator walking till the end of the strip that shall be collected
            auto xCollect2 = vMaxima.front().second;
            //save SoC index
            pSeeds->xStats.index_of_strip = uiSoCIndex++;
            // all these things are not used at the moment...
            pSeeds->xStats.uiInitialQueryBegin = xCollect2->start();
            pSeeds->xStats.uiInitialRefBegin = xCollect2->start_ref();
            pSeeds->xStats.uiInitialQueryEnd = xCollect2->end();
            pSeeds->xStats.uiInitialRefEnd = xCollect2->end_ref();
            nucSeqIndex end = getPositionForBucketing(uiQLen, *xCollect2) + uiStripSize;
            while(
                xCollect2 != vSeeds.end() &&
                end >= getPositionForBucketing(uiQLen, *xCollect2))
            {
                // save the beginning and end of the SoC
                // all these things are not used at the moment...
                if(xCollect2->start() < pSeeds->xStats.uiInitialQueryBegin)
                    pSeeds->xStats.uiInitialQueryBegin = xCollect2->start();
                if(xCollect2->start_ref() < pSeeds->xStats.uiInitialRefBegin)
                    pSeeds->xStats.uiInitialRefBegin = xCollect2->start_ref();
                if(xCollect2->end() > pSeeds->xStats.uiInitialQueryEnd)
                    pSeeds->xStats.uiInitialQueryEnd = xCollect2->end();
                if(xCollect2->end_ref() > pSeeds->xStats.uiInitialRefEnd)
                    pSeeds->xStats.uiInitialRefEnd = xCollect2->end_ref();
                pSeeds->xStats.num_seeds_in_strip++;
                assert(xCollect2->start() <= xCollect2->end());
                assert(xCollect2->end() <= pQuerySeq->length());
                //if the iterator is still within the strip add the seed and increment the iterator
                pSeeds->push_back(*(xCollect2++));
            }//while
            //move to the next strip
            std::pop_heap (vMaxima.begin(), vMaxima.end(), vHeapOrder); vMaxima.pop_back();

            return pSeeds;
        }// function

        inline void push_back_no_overlap(
                const SoCOrder& rCurrScore, 
                const std::vector<Seed>::iterator itStrip,
                const nucSeqIndex uiCurrStart,
                const nucSeqIndex uiCurrEnd,
            )
        {
            DEBUG(assert(!bInPriorityMode);)
            if(
                vMaxima.empty() || 
                uiLastEnd < uiCurrStart ||
                vMaxima.back().first < rCurrScore)
            {
                // if we reach this point we want to save the current SoC
                if( !vMaxima.empty() && uiLastEnd >= uiCurrStart)
                    // the new and the last SoC overlap and the new one has a higher score
                    // so we want to replace the last SoC
                    vMaxima.pop_back();
                // else the new and the last SoC do not overlapp
                // so we want to save the new SOC
                vMaxima.push_back(std::make_pair(rCurrScore, itStrip));
                // we need to remember the end position of the new SoC
                uiLastEnd = uiCurrEnd;
            }// if
            // else new and last SoC overlap and the new one has a lower score => ignore the new one
        }// function

        inline void make_heap()
        {
            DEBUG(assert(!bInPriorityMode);bInPriorityMode = true;)
            // make a max heap from the SOC starting points according to the scores, 
            // so that we can extract the best SOC first
            const auto vHeapOrder = 
                []
                (
                    const std::pair<SoCOrder, std::vector<Seed>::iterator>& rA,
                    const std::pair<SoCOrder, std::vector<Seed>::iterator>& rB
                )
                {
                    return rA.first < rB.first;
                }//lambda
            ;
            std::make_heap(vMaxima.begin(), vMaxima.end(), vHeapOrder);
        }// function
    }; //class
}// namspace libMA


void exportSoC();