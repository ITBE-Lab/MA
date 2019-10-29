/**
 * @file minimizerSeeding.cpp
 * @author Markus Schmidt
 */
#include "module/minimizerSeeding.h"

#ifdef WITH_ZLIB
using namespace libMA;

#ifdef WITH_PYTHON

void exportMinimizerSeeding( py::module& rxPyModuleId )
{
    // export the BinarySeeding class
    exportModule<MinimizerSeeding>( rxPyModuleId, "MinimizerSeeding" );
} // function
#endif
#endif