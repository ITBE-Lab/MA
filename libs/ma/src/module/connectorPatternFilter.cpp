/**
 * @file connectorPatternFilter.cpp
 * @author Markus Schmidt
 */
#include "module/connectorPatternFilter.h"

using namespace libMA;


#ifdef WITH_PYTHON

void exportConnectorPatternFilter( py::module& rxPyModuleId )
{
    // export the ConnectorPatternFilter class
    exportModule<ConnectorPatternFilter, std::shared_ptr<SV_DB>>( rxPyModuleId, "ConnectorPatternFilter" );
} // function
#endif