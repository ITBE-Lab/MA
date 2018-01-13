/** 
 * @file stripOfConsideration.h
 * @brief Implements a Strip of Consideration.
 * @author Markus Schmidt
 */
#ifndef STRIP_OF_CONSIDERATION_H
#define STRIP_OF_CONSIDERATION_H

#include "container/segment.h"
#include "module/module.h"
#include "module/needlemanWunsch.h"
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
        /// @brief Maximum ambiguity for a seed to be considered.
        unsigned int uiMaxAmbiguity = 500;
        /// @brief Minimum amount of seeds for a strip to be considered
        unsigned int minSeeds = 3;
        /// #brief the amount of SOCs to create
        unsigned int numStrips = 10;
        /// @brief Minimum nucleotides covered by seeds for a strip to be considered
        float minSeedLength = 0.1;
        /// @brief maximum amount of seeds total
        float dMaxSeeds = 7.f;

        /**
        * @brief skip seeds with too much ambiguity
        * @details
        * True: skip all seeds with to much ambiguity
        * False: use max_hits instances of the seeds with more ambiguity
        */
        bool bSkipLongBWTIntervals = true;
        
        static nucSeqIndex getPositionForBucketing(nucSeqIndex uiQueryLength, const Seed xS);

        static nucSeqIndex getStripSize(nucSeqIndex uiQueryLength);

        void forEachNonBridgingSeed(
                std::shared_ptr<SegmentVector> pVector,
                std::shared_ptr<FMIndex> pxFM_index,std::shared_ptr<Pack> pxRefSequence,
                std::shared_ptr<NucSeq> pxQuerySeq,
                std::function<void(Seed)> fDo,
                nucSeqIndex addSize// = 0 (default)
            );

        void sort(std::vector<Seed>& vSeeds, nucSeqIndex qLen);

    public:

        StripOfConsideration(){}//default constructor

        StripOfConsideration( 
                unsigned int uiMaxAmbiguity,
                unsigned int minSeeds,
                unsigned int numStrips,
                float minSeedLength,
                float dMaxSeeds
            )
                :
            uiMaxAmbiguity(uiMaxAmbiguity),
            minSeeds(minSeeds),
            numStrips(numStrips),
            minSeedLength(minSeedLength),
            dMaxSeeds(dMaxSeeds)
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