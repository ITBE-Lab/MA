/**
 * @file default_parameters.h
 * @brief Contains the default parameters for MA.
 * @author Markus Schmidt
 * @details
 * @todo this file should be removed
 */
#ifndef DEFAULT_PARAMETERS_H
#define DEFAULT_PARAMETERS_H

#include "support.h"
#include <string>

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
extern EXPORTED double dUnpaired; // @todo this should be a called paired bonus dPairedBonus
extern EXPORTED size_t uiMean; // @todo this should be a double -> dMean
extern EXPORTED double fStd;
extern EXPORTED size_t uiZDrop;
extern EXPORTED size_t uiSVPenalty;
extern EXPORTED int iBandwidthDPExtension;

} // namespace defaults
} // namespace libMA

#endif