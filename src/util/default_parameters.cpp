#include "util/default_parameters.h"

namespace libMA
{
    namespace defaults
    {        
        int iMatch = 2;      // score for a DP match (used in SoC width computation)
        int iMissMatch = 4;  // penalty for a DP missmatch
        int iGap = 4;       // penalty for a DP gap opening (used in SoC width computation)
        int iExtend = 2;     // penalty for a DP gap extension (used in SoC width computation)
        int iGap2 = 24;       // penalty for a DP gap opening (used in SoC width computation)
        int iExtend2 = 1;     // penalty for a DP gap extension (used in SoC width computation)
        size_t uiUnpaired = 17;  // penalty for unpaired reads
        size_t uiMean = 400;     // mean distance for paired reads
        double fStd = 150;      // standard deviation for distance of paired reads
        size_t uiReportN = 0;    // report n alignments
        size_t uiMaxAmbiguity = 100;     // maximal ambiguity of seeds
        size_t uiMinLen = 16;    // minimal seed length
        size_t uiMinAmbiguity = 0;       // stop the extension process if seeds are less ambiguous
        // @todo we effectively disabled this parameter
        size_t uiMinSeedSizeDrop = 15;   // minimum length for seeds to count towards the drop of
        size_t uiMaxTries = 30;          // maximal number of SoCs
        size_t uiMinTries = 2;          // minimal number of SoCs
        size_t uiMaxEqualScoreLookahead = 3;     // lookahead distance for short queries
        size_t uiSwitchQLen = 800;       // q len to switch between break criteria
        uint64_t uiMaxGapArea = 10000;     // break alignments in harmonization if gap is larger
        uint64_t uiPadding = 1000;          // padding for DP
        size_t uiSoCWidth = 0;          // set a fixed SoC width; 0 = use the formula
        bool bFindMode = false;        // true: don't do DP
        bool bOptimisticGapEstimation = true;  // how to estimate gap costs in harmonization
        bool bSkipLongBWTIntervals = true;     // pick samples from long SAintervals or skip them
        bool bNormalDist = false;     // use normal distribution to model gaps between paired reads
        bool bUniformDist = false;     // use normal distribution to model gaps between paired reads
        float fGiveUp = 0.002;  // do not store socs with score less than this * query len
        float fRelMinSeedSizeAmount = 0.005;    // minimum seed coverage to consider a query
        float fScoreDiffTolerance = 0.0001;     // break if the harm score falls faster than this
        float fSoCScoreMinimum = 0;             // minimum score used for the SoC width computation
        float fMinimalQueryCoverage = 1.1;      // does nothing at the moment @todo
        float fScoreTolerace = 0.1;             // break if the SoC score drops faster than this
        size_t uiCurrHarmScoreMin = 18;          // minimal score after the harmonization
        std::string sParameterSet = "fast";           // name of the used presetting (if any)
        std::string sSeedSet;             // name of the seed set that shall be computed
        // disable fGiveUp and fRelMinSeedSizeAmount for short genomes
        // @todo apply this parameter
        size_t uiGenomeSizeDisable = 10000000;
        bool bDisableHeuristics = false; // disable all heuristics in the harmonization
        // only output secondary alignments with a score larger than this * score of prim. alignment
        float fMinSecScoreRatio = .25;
        // can only be accessed by pacBio presetting at the moment...
        bool bDisableGapCostEstimationCutting = false; // do not remove seeds that are far away to generate alignments with positive score

        //artifact filter:
        // filter seeds if the difference between the delta distance to it's predecessor and
        // successor is less then <num> percent (set to 1 to disable filter)...
        double dMaxDeltaDist = 0.1;
        // AND the delta distance to it's pre- and successor is more than <num> nt.
        uint64_t uiMinDeltaDist = 16;
        // maximal amount of nucleotides primary and supplementary alignments can overlap.
        double dMaxOverlapSupplementary = 0.1;
        // each primary alignment can have at most x supplementary ones.
        size_t uiMaxSupplementaryPerPrim = 1;

        size_t uiSVPenalty = 100;
        
        // @todo
        double dMaxSVRatio = 0.01;
        int64_t iMinSVDistance = 500;
        
        size_t uiZDrop = 200;
        
#ifdef WITH_PYTHON
        void exportDefaults()
        {
            boost::python::def("configureFast", &configureFast);
            boost::python::def("configureAccurate", &configureAccurate);
        }//function
#endif
    };
};