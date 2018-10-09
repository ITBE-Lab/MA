#include "util/default_parameters.h"

namespace libMA
{
namespace defaults
{
int EXPORTED iMatch = 2; // score for a DP match (used in SoC width computation)
int EXPORTED iMissMatch = 4; // penalty for a DP missmatch
int EXPORTED iGap = 4; // penalty for a DP gap opening (used in SoC width computation)
int EXPORTED iExtend = 2; // penalty for a DP gap extension (used in SoC width computation)
int EXPORTED iGap2 = 24; // penalty for a DP gap opening (used in SoC width computation)
int EXPORTED iExtend2 = 1; // penalty for a DP gap extension (used in SoC width computation)
size_t EXPORTED uiUnpaired = 17; // penalty for unpaired reads
size_t EXPORTED uiMean = 400; // mean distance for paired reads
double EXPORTED fStd = 150; // standard deviation for distance of paired reads
size_t EXPORTED uiReportN = 0; // report n alignments
size_t EXPORTED uiMaxAmbiguity = 100; // maximal ambiguity of seeds
size_t EXPORTED uiMinLen = 16; // minimal seed length
size_t EXPORTED uiMinAmbiguity = 0; // stop the extension process if seeds are less ambiguous
// @todo we effectively disabled this parameter
size_t EXPORTED uiMinSeedSizeDrop = 15; // minimum length for seeds to count towards the drop of
size_t EXPORTED uiMaxTries = 30; // maximal number of SoCs
size_t EXPORTED uiMinTries = 1;//2; // minimal number of SoCs @todo
size_t EXPORTED uiMaxEqualScoreLookahead = 3; // lookahead distance for short queries
size_t EXPORTED uiSwitchQLen = 800; // q len to switch between break criteria
uint64_t EXPORTED uiMaxGapArea = 10000; // break alignments in harmonization if gap is larger
uint64_t EXPORTED uiPadding = 1000; // padding for DP
size_t EXPORTED uiSoCWidth = 0; // set a fixed SoC width; 0 = use the formula
bool EXPORTED bOptimisticGapEstimation = true; // how to estimate gap costs in harmonization
bool EXPORTED bSkipLongBWTIntervals = true; // pick samples from long SAintervals or skip them
bool EXPORTED bNormalDist = false; // use normal distribution to model gaps between paired reads
bool EXPORTED bUniformDist = false; // use normal distribution to model gaps between paired reads
float EXPORTED fGiveUp = 0.002f; // do not store socs with score less than this * query len
float EXPORTED fRelMinSeedSizeAmount = 0.005f; // minimum seed coverage to consider a query
float EXPORTED fScoreDiffTolerance = 0.0001f; // break if the harm score falls faster than this
float EXPORTED fSoCScoreMinimum = 0.0f; // minimum score used for the SoC width computation
float EXPORTED fMinimalQueryCoverage = 1.1f; // does nothing at the moment @todo
float EXPORTED fScoreTolerace = 0.1f; // break if the SoC score drops faster than this
size_t EXPORTED uiCurrHarmScoreMin = 18; // minimal score after the harmonization
std::string EXPORTED sParameterSet = "fast"; // name of the used presetting (if any)
std::string EXPORTED sSeedSet; // name of the seed set that shall be computed
// disable fGiveUp and fRelMinSeedSizeAmount for short genomes
size_t EXPORTED uiGenomeSizeDisable = 10000000;
bool EXPORTED bDisableHeuristics = false; // disable all heuristics in the harmonization
// only output secondary alignments with a score larger than this * score of prim. alignment
float EXPORTED fMinSecScoreRatio = .25;
// can only be accessed by pacBio presetting at the moment...
bool EXPORTED bDisableGapCostEstimationCutting =
    false; // do not remove seeds that are far away to generate alignments with positive score

// artifact filter:
// filter seeds if the difference between the delta distance to it's predecessor and
// successor is less then <num> percent (set to 1 to disable filter)...
double EXPORTED dMaxDeltaDist = 0.1;
// AND the delta distance to it's pre- and successor is more than <num> nt.
uint64_t EXPORTED uiMinDeltaDist = 16;
// maximal amount of nucleotides primary and supplementary alignments can overlap.
double EXPORTED dMaxOverlapSupplementary = 0.1;
// each primary alignment can have at most x supplementary ones.
size_t EXPORTED uiMaxSupplementaryPerPrim = 1;

size_t EXPORTED uiSVPenalty = 100;

// Minimal bandwith when filling the gap between seeds.
// This bandwidth gets increased if the seeds are not on a diagonal.
int EXPORTED iMinBandwidthGapFilling = 20;
// When extending the end of a alignment we use this bandwidth.
int EXPORTED iBandwidthDPExtension = 512;

// @todo
double EXPORTED dMaxSVRatio = 0.01;
int64_t EXPORTED iMinSVDistance = 500;

size_t EXPORTED uiZDrop = 200;

#ifdef WITH_PYTHON
void exportDefaults( )
{
    boost::python::def( "configureFast", &configureFast );
    boost::python::def( "configureAccurate", &configureAccurate );
} // function
#endif
}; // namespace defaults
}; // namespace libMA