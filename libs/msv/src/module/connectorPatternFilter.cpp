/**
 * @file connectorPatternFilter.cpp
 * @author Markus Schmidt
 */
#include "module/connectorPatternFilter.h"

using namespace libMSV;

#ifdef WITH_PYTHON

#include "container/sv_db/py_db_conf.h"

void exportConnectorPatternFilter( py::module& rxPyModuleId )
{
    // export the ConnectorPatternFilter class
    exportModule<ConnectorPatternFilter<DBCon>>( rxPyModuleId, "ConnectorPatternFilter" );

} // function
#endif