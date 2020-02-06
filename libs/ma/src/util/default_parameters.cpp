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
double EXPORTED dUnpaired = 1.25; // penalty for unpaired reads
size_t EXPORTED uiMean = 400; // mean distance for paired reads
double EXPORTED fStd = 150; // standard deviation for distance of paired reads
size_t EXPORTED uiSVPenalty = 100;
// When extending the end of a alignment we use this bandwidth.
int EXPORTED iBandwidthDPExtension = 512;
size_t EXPORTED uiZDrop = 200;

} // namespace defaults
} // namespace libMA