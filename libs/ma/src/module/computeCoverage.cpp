/**
 * @file computeCoverage.cpp
 * @author Markus Schmidt
 */
#include "module/computeCoverage.h"

using namespace libMA;


#ifdef WITH_PYTHON

void exportComputeCoverage( py::module& rxPyModuleId )
{
    // export the ComputeCoverage class
    exportModule<ComputeCoverage>( rxPyModuleId, "ComputeCoverage" );
} // function
#endif