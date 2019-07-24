/**
 * @file sweepSvJumps.cpp
 * @author Markus Schmidt
 */
#include "module/sweepSvJumps.h"

using namespace libMA;


#ifdef WITH_PYTHON
void exportSweepSvJump( py::module& rxPyModuleId )
{
    exportModule<CompleteBipartiteSubgraphSweep, std::shared_ptr<SV_DB>, std::shared_ptr<Pack>, int64_t>(
        rxPyModuleId, "CompleteBipartiteSubgraphSweep" );
} // function
#endif