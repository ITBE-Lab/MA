/**
 * @file export.cpp
 * @author Markus Schmidt
 */
#include "msv/container/svJump.h"
#include "msv/module/sweepSvJumps.h"
#include "msv/module/svJumpsFromSeeds.h"
#include "msv/module/connectorPatternFilter.h"
#include "msv/container/sv_db/svSchema.h"
#include "msv/module/count_k_mers.h"
#include "msv/util/rendererSpeedup.h"

using namespace libMA;


#ifdef WITH_PYTHON


PYBIND11_MODULE( libMSV, libMsvModule )
{
    py::module::import("libMS");
    libMS::SubmoduleOrganizer xOrganizer( libMsvModule );
    exportSVJump( xOrganizer );
    exportSvJumpsFromSeeds( xOrganizer );
    exportSweepSvJump( xOrganizer );
    exportConnectorPatternFilter( xOrganizer );
    exportSoCDbWriter( xOrganizer );
    exportCountKMers( xOrganizer );
    exportRendererSpeedUp(xOrganizer);
} // function

#endif
