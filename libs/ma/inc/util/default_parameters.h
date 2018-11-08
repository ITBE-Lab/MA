/**
 * @file default_parameters.h
 * @brief Contains the default parameters for MA.
 * @author Markus Schmidt
 */
#ifndef DEFAULT_PARAMETERS_H
#define DEFAULT_PARAMETERS_H

#include "util/support.h"
#include <string>

#ifdef WITH_PYTHON
    #ifdef BOOST_PYTHON
        #include <boost/python.hpp>
    #else
        #include <pybind11/stl_bind.h>
        namespace py = pybind11;
    #endif
#endif

namespace libMA
{
namespace defaults
{

extern EXPORTED int iMatch;
extern EXPORTED int iMissMatch;
extern EXPORTED int iGap;
extern EXPORTED int iExtend;
extern EXPORTED int iGap2;
extern EXPORTED int iExtend2;
extern EXPORTED size_t uiUnpaired;
extern EXPORTED size_t uiMean;
extern EXPORTED double fStd;
extern EXPORTED size_t uiReportN;
extern EXPORTED size_t uiMaxAmbiguity;
extern EXPORTED size_t uiMinLen;
extern EXPORTED size_t uiMinAmbiguity;
extern EXPORTED size_t uiMinSeedSizeDrop;
extern EXPORTED size_t uiMaxTries;
extern EXPORTED size_t uiMinTries;
extern EXPORTED size_t uiMaxEqualScoreLookahead;
extern EXPORTED size_t uiSwitchQLen;
extern EXPORTED uint64_t uiMaxGapArea;
extern EXPORTED uint64_t uiPadding;
extern EXPORTED size_t uiSoCWidth;
extern EXPORTED bool bFindMode;
extern EXPORTED bool bOptimisticGapEstimation;
extern EXPORTED bool bSkipLongBWTIntervals;
extern EXPORTED bool bNormalDist;
extern EXPORTED bool bUniformDist;
extern EXPORTED float fGiveUp;
extern EXPORTED float fRelMinSeedSizeAmount;
extern EXPORTED float fScoreDiffTolerance;
extern EXPORTED float fSoCScoreMinimum;
extern EXPORTED float fMinimalQueryCoverage;
extern EXPORTED float fScoreTolerace;
extern EXPORTED size_t uiCurrHarmScoreMin;
extern EXPORTED std::string sParameterSet;
extern EXPORTED std::string sSeedSet;
extern EXPORTED size_t uiGenomeSizeDisable;
extern EXPORTED bool bDisableHeuristics;
extern EXPORTED float fMinSecScoreRatio;
extern EXPORTED double dMaxDeltaDist;
extern EXPORTED uint64_t uiMinDeltaDist;
extern EXPORTED double dMaxOverlapSupplementary;
extern EXPORTED size_t uiMaxSupplementaryPerPrim;
extern EXPORTED bool bDisableGapCostEstimationCutting;
extern EXPORTED double dMaxSVRatio;
extern EXPORTED int64_t iMinSVDistance;
extern EXPORTED size_t uiZDrop;
extern EXPORTED size_t uiSVPenalty;
extern EXPORTED int iMinBandwidthGapFilling;
extern EXPORTED int iBandwidthDPExtension;

inline void configureAccurate( )
{
    sParameterSet = "acc";
    sSeedSet = "SMEMs";
} // function

inline void configureFast( )
{
    sParameterSet = "fast";
    sSeedSet = "maxSpan";
} // function

// inline void configurePacBio()
//{
//    sParameterSet = "pacBio";
//    sSeedSet = "maxSpan";
//    bDisableHeuristics = true;
//    uiReportN = 3;
//}// function

#ifdef WITH_PYTHON
#ifdef WITH_BOOST
void exportDefaults( );
#else
void exportDefaults( py::module& rxPyModuleId );
#endif
#endif
} // namespace defaults
} // namespace libMA

#endif
