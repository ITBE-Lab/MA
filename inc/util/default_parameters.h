/** 
 * @file default_parameters.h
 * @brief Contains the default parameters for MA.
 * @author Markus Schmidt
 */
#ifndef DEFAULT_PARAMETERS_H
#define DEFAULT_PARAMETERS_H

#include <string>

namespace libMA
{
    namespace defaults
    {
        std::string
            uiMatch = "3",
            uiMissMatch = "4",
            uiOpen = "6",
            uiExtend = "1",
            uiUnpaired = "17",
            uiMean = "400",
            uiStd = "150",
            uiReportN = "0",
            uiMaxAmbiguity,
            uiMinLen = "16",
            uiMinAmbiguity = "0",
            uiMinSeedSizeDrop = "15",
            uiMaxTries = "50",
            uiMaxEqualScoreLookahead = "3",
            uiSwitchQLen = "800",
            uiMaxGapArea = "10000",
            uiPadding = "500",
            bFindMode = "false",
            bOptimisticGapEstimation = "true",
            bSkipLongBWTIntervals = "true",
            fGiveUp = "0.002",
            fRelMinSeedSizeAmount = "0.005",
            fScoreDiffTolerance = "0.0001",
            fSoCScoreMinimum = "0",
            fMinimalQueryCoverage = "1.1",
            fScoreTolerace = "0.1",
            uiCurrHarmScoreMin = "18",
            sParameterSet,
            sSeedSet
        ;

        inline void configureAccurate()
        {
            sParameterSet = "accurate";
            sSeedSet = "SMEMs";
            uiMaxAmbiguity = "100";
        }// function
        inline void configureFast()
        {
            sParameterSet = "fast";
            sSeedSet = "maxSpanning";
            uiMaxAmbiguity = "1000";
        }// function
    }// namespace
}// namespace

#endif
