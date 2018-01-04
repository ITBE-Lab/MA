/** 
 * @file stripOfConsideration.h
 * @brief Implements a Strip of Consideration.
 * @author Markus Schmidt
 */
#ifndef STRIP_OF_CONSIDERATION_H
#define STRIP_OF_CONSIDERATION_H

#include "container/segment.h"
#include "module/module.h"
#include <cmath>

namespace libMABS
{

    /**
     * @brief Used to quickly find areas with high density of @ref Seed "seeds".
     * @ingroup module
     */
    class StripOfConsideration: public Module
    {
    public:
        /// @brief The strip of consideration size.
        nucSeqIndex uiStripSize = 10000;
        /// @brief Maximum ambiguity for a seed to be considered.
        unsigned int uiMaxAmbiguity = 500;
        /// @brief Minimum amount of seeds for a strip to be considered
        unsigned int minSeeds = 3;
        /// @brief Minimum nucleotides covered by seeds for a strip to be considered
        float minSeedLength = 0.1;
        /// @brief maximum amount of seeds total
        float dMaxSeeds = 7.f;
        /// @brief maximum amount of seeds total
        float dMaxSeeds2 = 7.f;

        /**
        * @brief skip seeds with too much ambiguity
        * @details
        * True: skip all seeds with to much ambiguity
        * False: use max_hits instances of the seeds with more ambiguity
        */
        bool bSkipLongBWTIntervals = true;
        
    private:
        inline nucSeqIndex getPositionForBucketing(nucSeqIndex uiQueryLength, const Seed xS) const 
        { 
            return xS.start_ref() + (uiQueryLength - xS.start()); 
        }//function

        void forEachNonBridgingSeed(
                std::shared_ptr<SegmentVector> pVector,
                std::shared_ptr<FMIndex> pxFM_index,std::shared_ptr<Pack> pxRefSequence,
                std::shared_ptr<NucSeq> pxQuerySeq,
                std::function<void(Seed)> fDo,
                nucSeqIndex addSize// = 0 (default)
            );

        void linearSort(std::vector<std::tuple<Seed, bool>>& vSeeds, nucSeqIndex qLen)
        {
            unsigned int max_bits_used = 34;
            unsigned int n = vSeeds.size();
    
            if( 2*max_bits_used*n / log2(n) < n*log2(n) )
            {
                //quicksort is faster
                std::sort(
                    vSeeds.begin(), vSeeds.end(),
                    [&]
                    (const std::tuple<Seed, bool> a, const std::tuple<Seed, bool> b)
                    {
                        return getPositionForBucketing(qLen, std::get<0>(a)) 
                                < getPositionForBucketing(qLen, std::get<0>(b));
                    }//lambda
                );//sort function call
                return;
            }//if
            //radix sort is faster than quicksort

            unsigned int amount_buckets = max_bits_used/log2(n);
            unsigned int iter = 0;
            std::vector<std::vector<std::tuple<Seed, bool>>> xBuckets1(amount_buckets);
            std::vector<std::vector<std::tuple<Seed, bool>>> xBuckets2(amount_buckets);
            auto pBucketsNow = &xBuckets1;
            auto pBucketsLast = &xBuckets2;
            for(auto& xSeed : vSeeds)
                (*pBucketsLast)[0].push_back(xSeed);

            //2**max_bits_used maximal possible genome size
            while(std::pow(amount_buckets,iter) <= std::pow(2,max_bits_used))
            {
                for(auto& xBucket : *pBucketsNow)
                    xBucket.clear();
                for(auto& xBucket : *pBucketsLast)
                    for(auto& xSeed : xBucket)
                    {
                        nucSeqIndex index = getPositionForBucketing(qLen, std::get<0>(xSeed));
                        index = (nucSeqIndex)( index / std::pow(amount_buckets,iter) ) % amount_buckets;
                        (*pBucketsNow)[index].push_back(xSeed);
                    }//for
                iter++;
                //swap last and now
                auto pBucketsTemp = pBucketsNow;
                pBucketsNow = pBucketsLast;
                pBucketsLast = pBucketsTemp;
            }//while
            vSeeds.clear();
            for(auto& xBucket : *pBucketsLast)
                for(auto& xSeed : xBucket)
                    vSeeds.push_back(xSeed);
        }//function

    public:

        StripOfConsideration(){}//default constructor

        StripOfConsideration(
                nucSeqIndex uiStripSize, 
                unsigned int uiMaxAmbiguity,
                unsigned int minSeeds,
                float minSeedLength,
                float dMaxSeeds,
                float dMaxSeeds2
            )
                :
            uiStripSize(uiStripSize),
            uiMaxAmbiguity(uiMaxAmbiguity),
            minSeeds(minSeeds),
            minSeedLength(minSeedLength),
            dMaxSeeds(dMaxSeeds),
            dMaxSeeds2(dMaxSeeds2)
        {}//constructor

        std::shared_ptr<Container> execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - SegmentVector
         * - Seeds
         * - NucSeq
         * - Pack
         * - FMIndex
         */
        ContainerVector getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - ContainerVector(Seeds)
         */
        std::shared_ptr<Container> getOutputType() const;

        std::string getName() const
        {
            return "StripOfConsideration";
        }
    };//class
}//namspace libMABS

/**
 * @brief export the bucketing @ref Module "module" to python.
 * @ingroup export
 */
void exportStripOfConsideration();

#endif