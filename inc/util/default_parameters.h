/** 
 * @file default_parameters.h
 * @brief Contains the default parameters for MA.
 * @author Markus Schmidt
 */
#ifndef DEFAULT_PARAMETERS_H
#define DEFAULT_PARAMETERS_H

#include <string>

#ifdef WITH_PYTHON
#include <boost/python.hpp>
#endif

namespace libMA
{
    namespace defaults
    {
        extern int iMatch;
        extern int iMissMatch;
        extern int iGap;
        extern int iExtend;
        extern int iGap2;
        extern int iExtend2;
        extern size_t uiUnpaired;
        extern size_t uiMean;
        extern double fStd;
        extern size_t uiReportN;
        extern size_t uiMaxAmbiguity;
        extern size_t uiMinLen;
        extern size_t uiMinAmbiguity;
        extern size_t uiMinSeedSizeDrop;
        extern size_t uiMaxTries;
        extern size_t uiMinTries;
        extern size_t uiMaxEqualScoreLookahead;
        extern size_t uiSwitchQLen;
        extern uint64_t uiMaxGapArea;
        extern uint64_t uiPadding;
        extern size_t uiSoCWidth;
        extern bool bFindMode;
        extern bool bOptimisticGapEstimation;
        extern bool bSkipLongBWTIntervals;
        extern bool bNormalDist;
        extern bool bUniformDist;
        extern float fGiveUp;
        extern float fRelMinSeedSizeAmount;
        extern float fScoreDiffTolerance;
        extern float fSoCScoreMinimum;
        extern float fMinimalQueryCoverage;
        extern float fScoreTolerace;
        extern size_t uiCurrHarmScoreMin;
        extern std::string sParameterSet;
        extern std::string sSeedSet;
        extern size_t uiGenomeSizeDisable;
        extern bool bDisableHeuristics;
        extern float fMinSecScoreRatio;
        extern double dMaxDeltaDist;
        extern uint64_t uiMinDeltaDist;
        extern double dMaxOverlapSupplementary;
        extern size_t uiMaxSupplementaryPerPrim;
        extern bool bDisableGapCostEstimationCutting;
        extern double dMaxSVRatio;
        extern int64_t iMinSVDistance;
        extern size_t uiZDrop;
        extern size_t uiSVPenalty;

        inline void configureAccurate()
        {
            sParameterSet = "acc";
            sSeedSet = "SMEMs";
        }// function

        inline void configureFast()
        {
            sParameterSet = "fast";
            sSeedSet = "maxSpan";
        }// function

        //inline void configurePacBio()
        //{
        //    sParameterSet = "pacBio";
        //    sSeedSet = "maxSpan";
        //    bDisableHeuristics = true;
        //    uiReportN = 3;
        //}// function

#ifdef WITH_PYTHON
        void exportDefaults();
#endif
    }// namespace
}// namespace 

#endif
