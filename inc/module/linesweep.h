/** 
 * @file linesweep.h
 * @brief Implements the linesweep @ref Module "module"
 * @author Markus Schmidt
 */
#ifndef LINESWEEP_H
#define LINESWEEP_H

#include "container/segment.h"
#include "container/soc.h"
#include "module/module.h"

#define USE_RANSAC ( 1 )

namespace libMA
{
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
         * @brief The shadow of a Seed.
         * @details
         * Each perfect match "casts a shadow" at the left and right border of the strip.
         * Each shadow is stored in one of these data structures.
         */
        class ShadowInterval: public Interval<int64_t>{
        public:
            Seeds::iterator pSeed;

            /**
             * @brief Creates a new shadow.
             * @details
             * The linesweep algorithm disables seeds.
             * Therefore the iterator is required in order to delete the respective seed from its list.
             */
            ShadowInterval(
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
            ShadowInterval( const ShadowInterval& rOther )
                    :
                Interval(rOther),
                pSeed(rOther.pSeed)
            {}//copy constructor

            bool within(const ShadowInterval& rOther)
            {
                return start() >= rOther.start() && end() <= rOther.end();
            }//function
        };//class


        /**
        * @brief Implements the linesweep algorithm.
        * @details
        * The algorithm has to be run on left and right shadows,
        * therefore it is provided as individual function.
        */
        std::shared_ptr<std::vector<std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>>>
        EXPORTED linesweep(
            std::shared_ptr<std::vector<
                std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>
            >> pShadows,
            const int64_t uiRStart,
            const double fAngle
        );
        
        inline double deltaDistance(
                const Seed& rSeed,
                const double fAngle,
                const int64_t uiRStart
            )
        {
            #define PI 3.14159265

            double y = rSeed.start_ref() + rSeed.start() / std::tan( PI/2 - fAngle );
            double x = ( y - uiRStart ) * std::sin(fAngle);
            double x_1 = rSeed.start() / std::sin( PI/2 - fAngle);

            return std::abs(x - x_1);
        }// method

        std::shared_ptr<Seeds> applyFilters(std::shared_ptr<Seeds>& pIn) const;

    public:
        /**
         * @brief true: estimate all possible position as matches in the gap cost filter
         * @details
         * After the linesweep picks a subset of all seeds so that the overall score
         * becomes optimal.
         * This filter has a linear time complexity.
         * It estimates the penalty for gaps between seeds.
         * We have two options here:
         * True -> the gapcost is estimated optimistically (as small as possible)
         * False -> we assume that the score for matches/missmatches is roughly equal within the gap
         * 
         * @note not beeing optimistic here has negative affect on the accuracy
         * but improves runtime significantly
         */
        bool optimisticGapEstimation = defaults::bOptimisticGapEstimation;

        /// @brief If the seeds cover less that x percent of the query we use SW, 
        /// otherwise we fill in the gaps.
        double fMinimalQueryCoverage = defaults::fMinimalQueryCoverage;// 0.25 before 1.05.18 however python disables this

        /// @brief Stop the SoC extraction if the harmonized score drops more than x.
        double fScoreTolerace = defaults::fScoreTolerace;

        /// @brief Extract at most x SoCs.
        unsigned int uiMaxTries = defaults::bFindMode ? defaults::uiReportN : defaults::uiMaxTries;

        /// @brief Lookahead for the equality break criteria.
        unsigned int uiMaxEqualScoreLookahead = defaults::uiMaxEqualScoreLookahead;

        /// @brief Consider two scores equal if they do not differ by more than x (relative to the total score).
        float fScoreDiffTolerance = defaults::fScoreDiffTolerance;

        /// @brief switch between the two break criteria based on weather the query len is larger than x.
        nucSeqIndex uiSwitchQLen = defaults::uiSwitchQLen;

        /// @brief minimum accumulated seed length after the harmonization as absolute value.
        nucSeqIndex uiCurrHarmScoreMin = defaults::uiCurrHarmScoreMin;
        /// @brief minimum accumulated seed length after the harmonization relative to query length.
        float fCurrHarmScoreMinRel = defaults::fGiveUp;

        bool bDoHeuristics = !defaults::bDisableHeuristics;

        bool bDoGapCostEstimationCutting = !defaults::bDisableGapCostEstimationCutting;

        double dMaxDeltaDist = defaults::dMaxDeltaDist;

        nucSeqIndex uiMinDeltaDist = defaults::uiMinDeltaDist;

        double dMaxSVRatio = defaults::dMaxSVRatio;

        int64_t iMinSVDistance = defaults::iMinSVDistance;
        
        nucSeqIndex uiMaxGapArea = defaults::uiMaxGapArea;

        size_t uiSVPenalty = defaults::uiSVPenalty;

        size_t uiUONAccMaxSize = defaults::uiUONAccMaxSize;

        LinearLineSweep() {}//default constructor

        //overload
        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> pInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - Seeds
         */
        ContainerVector EXPORTED getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - Seeds
         */
        std::shared_ptr<Container> EXPORTED getOutputType() const;

        std::string getName() const
        {
            return "LineSweep2";
        }

        std::string getFullDesc() const
        {
            return std::string("LineSweep2(") + 
                std::to_string(optimisticGapEstimation) + "," +
                std::to_string(fMinimalQueryCoverage) + "," +
                std::to_string(fScoreTolerace) + "," +
                std::to_string(uiMaxEqualScoreLookahead) + "," +
                std::to_string(fScoreDiffTolerance) + "," +
                std::to_string(uiMaxTries) + ")"
                ;
        }//function
    };//class
}//namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief Exposes the LineSweep @ref Module "module" to boost python.
 * @ingroup export
 */
void exportLinesweep();
#endif

#endif