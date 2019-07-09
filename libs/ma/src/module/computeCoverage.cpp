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
    exportModule<ComputeCoverage,std::shared_ptr<SV_DB>, int64_t, int64_t>( rxPyModuleId, "ComputeCoverage" );
} // function
#endif