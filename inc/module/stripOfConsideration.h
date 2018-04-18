/** 
 * @file stripOfConsideration.h
 * @brief Implements a Strip of Consideration.
 * @author Markus Schmidt
 */
#ifndef STRIP_OF_CONSIDERATION_H
#define STRIP_OF_CONSIDERATION_H

#include "container/segment.h"
#include "container/soc.h"
#include "module/module.h"
#include <cmath>

namespace libMA
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
        /// @brief the amount of SOCs to create
        unsigned int numStrips = 10;
        /**
         * @brief Minimal SoC score.
         * @details
         * Must be [0,-inf]!
         * The minimal score that should be allowed during SoC collection.
         * Is interpreted relative to the query length.
         */
        const float fScoreMinimum = 0;
        /**
         * @brief If the best SoC has seeds of accumulative length smaller than this, abort.
         * @details
         * Is multiplied by query length.
         * 0 = never abort.
         */
        const float fGiveUp = 0;

        /**
         * @todo
         */
        const float fIncludeAmbiguity = 0;

        /**
        * @brief skip seeds with too much ambiguity
        * @details
        * True: skip all seeds with to much ambiguity
        * False: use max_hits instances of the seeds with more ambiguity
        */
        bool bSkipLongBWTIntervals = true;
        
        inline static nucSeqIndex getPositionForBucketing(nucSeqIndex uiQueryLength, const Seed& xS)
        { 
            return xS.start_ref() + (uiQueryLength - xS.start()); 
        }//function

        inline nucSeqIndex getStripSize(
                nucSeqIndex uiQueryLength,
                int iMatch,
                int iExtend,
                int iGap
            ) const
        {
            return 100;
            return (iMatch * uiQueryLength - iGap) / iExtend - (int64_t)(fScoreMinimum * uiQueryLength);
        }//function

        void EXPORTED forEachNonBridgingSeed(
                std::shared_ptr<SegmentVector> pVector,
                std::shared_ptr<FMIndex> pxFM_index,std::shared_ptr<Pack> pxRefSequence,
                std::shared_ptr<NucSeq> pxQuerySeq,
                std::function<void(Seed)> fDo,
                nucSeqIndex addSize// = 0 (default)
            );

    public:

        StripOfConsideration(){}//default constructor

        StripOfConsideration( 
                unsigned int uiMaxAmbiguity,
                unsigned int numStrips,
                float fScoreMinimum,
                float fGiveUp,
                float fIncludeAmbiguity
            )
                :
            uiMaxAmbiguity(uiMaxAmbiguity),
            numStrips(numStrips),
            fScoreMinimum(fScoreMinimum),
            fGiveUp(fGiveUp),
            fIncludeAmbiguity(fIncludeAmbiguity)
        {}//constructor

        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> vpInput);

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
        ContainerVector EXPORTED getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - ContainerVector(Seeds)
         */
        std::shared_ptr<Container> EXPORTED getOutputType() const;

        std::string getName() const
        {
            return "StripOfConsideration";
        }

        std::string getFullDesc() const
        {
            return "StripOfConsideration(" + 
                std::to_string(uiMaxAmbiguity) + "," +
                std::to_string(fScoreMinimum) + "," +
                std::to_string(fGiveUp) + "," +
                std::to_string(fIncludeAmbiguity) + "," +
                std::to_string(numStrips) + ")";
        }//function
    };//class

    /**
     * @brief Used to quickly find areas with high density of @ref Seed "seeds".
     * @ingroup module
     */
    class StripOfConsideration2: public Module
    {
    public:
        /// @brief Maximum ambiguity for a seed to be considered.
        unsigned int uiMaxAmbiguity = 500;
        /// @brief the amount of SOCs to create
        unsigned int numStrips = 10;

        static nucSeqIndex EXPORTED getStripSize(nucSeqIndex uiQueryLength);

        static nucSeqIndex EXPORTED getPositionForBucketing(nucSeqIndex uiQueryLength, const Seed xS);

    public:

        StripOfConsideration2(){}//default constructor

        StripOfConsideration2( 
                unsigned int uiMaxAmbiguity,
                unsigned int numStrips
            )
                :
            uiMaxAmbiguity(uiMaxAmbiguity),
            numStrips(numStrips)
        {}//constructor

        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> vpInput);

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
        ContainerVector EXPORTED getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - ContainerVector(Seeds)
         */
        std::shared_ptr<Container> EXPORTED getOutputType() const;

        std::string getName() const
        {
            return "StripOfConsideration2";
        }
    };//class
}//namspace libMA

/**
 * @brief export the bucketing @ref Module "module" to python.
 * @ingroup export
 */
void exportStripOfConsideration();

#endif