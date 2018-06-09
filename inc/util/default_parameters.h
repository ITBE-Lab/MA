/** 
 * @file default_parameters.h
 * @brief Contains the default parameters for MA.
 * @author Markus Schmidt
 */
#ifndef DEFAULT_PARAMETERS_H
#define DEFAULT_PARAMETERS_H

#include <string>


#define DEFAULT 0
#define EXPERIMENTAL 1

#define PARAMETER_SET DEFAULT


namespace libMA
{
    namespace defaults
    {
#if(PARAMETER_SET == DEFAULT)
        //@todo comment this file
        std::string
            uiMatch = "3",      // score for a DP match (used in SoC width computation)
            uiMissMatch = "4",  // penalty for a DP missmatch
            uiOpen = "6",       // penalty for a DP gap opening (used in SoC width computation)
            uiExtend = "1",     // penalty for a DP gap extension (used in SoC width computation)
            uiUnpaired = "17",  // penalty for unpaired reads
            uiMean = "400",     // mean distance for paired reads
            uiStd = "150",      // standard deviation for distance of paired reads
            uiReportN = "0",    // report n alignments
            uiMaxAmbiguity = "100",     // maximal ambiguity of seeds
            uiMinLen = "16",    // minimal seed length
            uiMinAmbiguity = "0",       // stop the extension process if seeds are less ambiguous
            // @todo we effectively disabled this parameter
            uiMinSeedSizeDrop = "15",   // minimum length for seeds to count towards the drop of
            uiMaxTries = "50",          // maximal number of SoCs
            uiMaxEqualScoreLookahead = "3",     // lookahead distance for short queries
            uiSwitchQLen = "800",       // q len to switch between break criteria
            uiMaxGapArea = "10000",     // break alignments in harmonization if gap is larger
            uiPadding = "500",          // padding for DP
            bFindMode = "false",        // true: don't do DP
            bOptimisticGapEstimation = "true",  // how to estimate gap costs in harmonization
            bSkipLongBWTIntervals = "true",     // pick samples from long SAintervals or skip them
            fGiveUp = "0.002",  // do not store socs with score less than this * query len
            fRelMinSeedSizeAmount = "0.005",    // minimum seed coverage to consider a query
            fScoreDiffTolerance = "0.0001",     // break if the harm score falls faster than this
            fSoCScoreMinimum = "0",             // minimum score used for the SoC width computation
            fMinimalQueryCoverage = "1.1",      // does nothing at the moment @todo
            fScoreTolerace = "0.1",             // break if the SoC score drops faster than this
            uiCurrHarmScoreMin = "18",          // minimal score after the harmonization
            sParameterSet,                      // name of the used presetting (if any)
            sSeedSet,                           // name of the seed set that shall be computed
            // disable fGiveUp and fRelMinSeedSizeAmount for short genomes
            // @todo apply this parameter
            uiGenomeSizeDisable = "10000000"
        ;
#elif(PARAMETER_SET == EXPERIMENTAL)
        std::string
            uiMatch = "3",      // score for a DP match (used in SoC width computation)
            uiMissMatch = "4",  // penalty for a DP missmatch
            uiOpen = "6",       // penalty for a DP gap opening (used in SoC width computation)
            uiExtend = "1",     // penalty for a DP gap extension (used in SoC width computation)
            uiUnpaired = "17",  // penalty for unpaired reads
            uiMean = "400",     // mean distance for paired reads
            uiStd = "150",      // standard deviation for distance of paired reads
            uiReportN = "0",    // report n alignments
            uiMaxAmbiguity = "1000",     // maximal ambiguity of seeds
            uiMinLen = "0",    // minimal seed length
            uiMinAmbiguity = "0",       // stop the extension process if seeds are less ambiguous
            // @todo we effectively disabled this parameter
            uiMinSeedSizeDrop = "0",   // minimum length for seeds to count towards the drop of
            uiMaxTries = "500",          // maximal number of SoCs
            uiMaxEqualScoreLookahead = "0",     // lookahead distance for short queries
            uiSwitchQLen = "0",       // q len to switch between break criteria
            uiMaxGapArea = "10000",     // break alignments in harmonization if gap is larger
            uiPadding = "500",          // padding for DP
            bFindMode = "false",        // true: don't do DP
            bOptimisticGapEstimation = "true",  // how to estimate gap costs in harmonization
            bSkipLongBWTIntervals = "true",     // pick samples from long SAintervals or skip them
            fGiveUp = "0",  // do not store socs with score less than this * query len
            fRelMinSeedSizeAmount = "1.1",    // minimum seed coverage to consider a query
            fScoreDiffTolerance = "0",     // break if the harm score falls faster than this
            fSoCScoreMinimum = "0",             // minimum score used for the SoC width computation
            fMinimalQueryCoverage = "1.1",      // does nothing at the moment @todo
            fScoreTolerace = "0",             // break if the SoC score drops faster than this
            uiCurrHarmScoreMin = "0",          // minimal score after the harmonization
            sParameterSet,                      // name of the used presetting (if any)
            sSeedSet,                           // name of the seed set that shall be computed
            // disable fGiveUp and fRelMinSeedSizeAmount for short genomes
            // @todo apply this parameter
            uiGenomeSizeDisable = "10000000"
        ;
#endif

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
    }// namespace
}// namespace

#endif
