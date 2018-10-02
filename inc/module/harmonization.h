/**
 * @file harmonization.h
 * @brief Implements the harmonization @ref Module "module"
 * @author Markus Schmidt
 */
#ifndef LINESWEEP_H
#define LINESWEEP_H

#include "container/segment.h"
#include "container/soc.h"
#include "module/module.h"

#define USE_RANSAC ( 1 )
#define PI 3.14159265

namespace libMA
{
/**
 * @brief Implements the Harmonization algorithm.
 * @ingroup module
 * @details
 * Removes all contradicting seeds.
 * This should only be used in combination with the StripOfConsideration module.
 */
class Harmonization : public Module<ContainerVector<std::shared_ptr<Seeds>>, false, SoCPriorityQueue, NucSeq>
{
  private:
    /**
     * @brief The shadow of a Seed.
     * @details
     * Each perfect match "casts a shadow" at the left and right border of the strip.
     * Each shadow is stored in one of these data structures.
     */
    class ShadowInterval : public Interval<int64_t>
    {
      public:
        Seeds::iterator pSeed;

        /**
         * @brief Creates a new shadow.
         * @details
         * The linesweep algorithm disables seeds.
         * Therefore the iterator is required in order to delete the respective seed from its list.
         */
        ShadowInterval( int64_t iBegin, int64_t iSize, Seeds::iterator pSeed )
            : Interval( iBegin, iSize ), pSeed( pSeed )
        {} // constructor

        /**
         * @brief Copy constructor
         */
        ShadowInterval( const ShadowInterval& rOther ) : Interval( rOther ), pSeed( rOther.pSeed )
        {} // copy constructor

        bool within( const ShadowInterval& rOther )
        {
            return start( ) >= rOther.start( ) && end( ) <= rOther.end( );
        } // function
    }; // class


    /**
     * @brief Implements the linesweep algorithm.
     * @details
     * The algorithm has to be run on left and right shadows,
     * therefore it is provided as individual function.
     */
    std::shared_ptr<std::vector<std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>>> EXPORTED
    linesweep( std::shared_ptr<std::vector<std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>>> pShadows,
               const int64_t uiRStart, const double fAngle );

    inline double deltaDistance( const Seed& rSeed, const double fAngle, const int64_t uiRStart )
    {
        double y = rSeed.start_ref( ) + rSeed.start( ) / std::tan( PI / 2 - fAngle );
        double x = ( y - uiRStart ) * std::sin( fAngle );
        double x_1 = rSeed.start( ) / std::sin( PI / 2 - fAngle );

        return std::abs( x - x_1 );
    } // method

    std::shared_ptr<Seeds> applyFilters( std::shared_ptr<Seeds>& pIn ) const;

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
    double fMinimalQueryCoverage = defaults::fMinimalQueryCoverage;

    /// @brief Stop the SoC extraction if the harmonized score drops more than x.
    double fScoreTolerace = defaults::fScoreTolerace;

    /// @brief Extract at most x SoCs.
    unsigned int uiMaxTries = defaults::bFindMode ? defaults::uiReportN : defaults::uiMaxTries;
    /// @brief Extract at least x SoCs.
    size_t uiMinTries = defaults::uiMinTries;

    /// @brief Lookahead for the equality break criteria.
    unsigned int uiMaxEqualScoreLookahead = defaults::uiMaxEqualScoreLookahead;

    /// @brief Consider two scores equal if they do not differ by more than x (relative to the total
    /// score).
    float fScoreDiffTolerance = defaults::fScoreDiffTolerance;

    /// @brief switch between the two break criteria based on weather the query len is larger than
    /// x.
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

    Harmonization( )
    {} // default constructor

    // overload
    virtual std::shared_ptr<ContainerVector<std::shared_ptr<Seeds>>>
        EXPORTED execute( std::shared_ptr<SoCPriorityQueue> pSoCIn, std::shared_ptr<NucSeq> pQuery );

}; // class
} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief Exposes the Harmonization @ref Module "module" to boost python.
 * @ingroup export
 */
void exportHarmonization( );
#endif

#endif