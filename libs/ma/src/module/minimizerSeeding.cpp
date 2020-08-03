/**
 * @file minimizerSeeding.cpp
 * @author Markus Schmidt
 */
#include "ma/module/minimizerSeeding.h"

#ifdef WITH_ZLIB
using namespace libMA;
using namespace libMS;

#ifdef WITH_PYTHON

void exportMinimizerSeeding( libMS::SubmoduleOrganizer& xOrganizer )
{
    // export the MinimizerSeeding class
    exportModule<MinimizerSeeding>( xOrganizer, "MinimizerSeeding" );
} // function
#endif
#endif