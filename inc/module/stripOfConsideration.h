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
        unsigned int uiMaxAmbiguity = 1000;
        /// @brief Minimum seed length.
        unsigned int uiMinLen = 16;
        /**
         * @brief Minimal SoC score.
         * @details
         * Must be [0,-inf]!
         * The minimal score that should be allowed during SoC collection.
         * Is interpreted relative to the query length.
         */
        float fScoreMinimum = 0;
        /**
         * @brief If the best SoC has seeds of accumulative length smaller than this, abort.
         * @details
         * Is multiplied by query length.
         * 0 = never abort.
         */
        const float fGiveUp = 0;
        unsigned int uiCurrHarmScoreMin = 18;
        
        /**
         * @brief disable fGiveUp and fRelMinSeedSizeAmount if genome is too short
         */
        const nucSeqIndex uiMinGenomeSize = 0;

        /**
        * @brief skip seeds with too much ambiguity
        * @details
        * True: skip all seeds with to much ambiguity
        * False: use max_hits instances of the seeds with more ambiguity
        */
        bool bSkipLongBWTIntervals = true;

        inline static nucSeqIndex getPositionForBucketing( nucSeqIndex uiQueryLength, 
                                                           const Seed& xS )
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
            return (iMatch * uiQueryLength - iGap) / iExtend - (int64_t)(fScoreMinimum * uiQueryLength);
        }//function

        template <class FUNCTOR>
        inline void forEachNonBridgingSeed(
                SegmentVector &rSegmentVector,
                FMIndex &rxFM_index,
                Pack &rxRefSequence,
                FUNCTOR&& fDo, //std::function<void(const Seed &rxS)> fDo,
                nucSeqIndex addSize = 0
            )
        {
            rSegmentVector.forEachSeed(
                rxFM_index, uiMaxAmbiguity, bSkipLongBWTIntervals,
                [&rxRefSequence, &rxFM_index, &fDo, &addSize]
                (Seed &&xS)
                {
                    // check if the match is bridging the forward/reverse strand 
                    // or bridging between two chromosomes
                    if ( !rxRefSequence.bridgingSubsection(
                            //prevent negative index
                            xS.start_ref() > addSize ? xS.start_ref() - addSize : 0,//from
                            //prevent index larger than reference
                            xS.end_ref() + addSize <= rxFM_index.getRefSeqLength() ?
                                xS.size() - 1 + addSize :
                                rxFM_index.getRefSeqLength() - xS.start_ref() - 1//to
                            ) 
                        )
                    {
                        //if non-bridging use this seed
                        fDo( xS );
                    }//if
                    //returning true since we want to continue extracting seeds
                    return true;
                }//lambda
            );
        }//method
        

        inline void emplaceAllNonBridgingSeed(
                SegmentVector &rSegmentVector,
                FMIndex &rxFM_index,
                Pack &rxRefSequence,
                std::vector<Seed> &rvSeedVector,
                const nucSeqIndex uiQLen
            )
        {
            rSegmentVector.emplaceAllEachSeeds(
                rxFM_index, uiMaxAmbiguity, uiMinLen, rvSeedVector,
                [&rxRefSequence, &rxFM_index, &rvSeedVector, &uiQLen]
                ()
                {
                    /*
                     * @note this bridging check is not required since we check weather a SoC
                     * is brigding in general.
                     * If any of the seeds within a SoC are bridging then the SoC is bridging.
                     * Could turn this into a debug assertion...
                     */
# if 0 // enable / disable the bridging check for all seeds
                    constexpr const nucSeqIndex addSize = 0;
                    auto& rS = rvSeedVector.back();
                    // check if the match is bridging the forward/reverse strand 
                    // or bridging between two chromosomes
                    if ( rxRefSequence.bridgingSubsection(
                            //prevent negative index
                            rS.start_ref() > addSize ? rS.start_ref() - addSize : 0,//from
                            //prevent index larger than reference
                            rS.end_ref() + addSize <= rxFM_index.getRefSeqLength() ?
                                rS.size() - 1 + addSize :
                                rxFM_index.getRefSeqLength() - rS.start_ref() - 1//to
                            ) 
                        )
                    {
                        //if bridging remove this seed
                        rvSeedVector.pop_back();
                    }//if
#if DELTA_CACHE == ( 1 )
                    else
                    {
                        // set the delta cache
                        rS.uiDelta = getPositionForBucketing( uiQLen, rS );
                    }// else
#endif
#elif DELTA_CACHE == ( 1 )
                    rvSeedVector.back().uiDelta = 
                        getPositionForBucketing( uiQLen, rvSeedVector.back() );
#endif
                    //returning true since we want to continue extracting seeds
                    return true;
                }//lambda
            );
        }//method

    public:

        StripOfConsideration(){}//default constructor

        StripOfConsideration(
                float fGiveUp
            )
                :
            fGiveUp(fGiveUp)
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
                std::to_string(fGiveUp) + "," + ")";
        }//function
    };//class


}//namspace libMA

#ifdef WITH_PYTHON
/**
 * @brief export the bucketing @ref Module "module" to python.
 * @ingroup export
 */
void exportStripOfConsideration();
#endif

#endif