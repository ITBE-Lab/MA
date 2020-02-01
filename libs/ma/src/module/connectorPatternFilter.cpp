/**
 * @file connectorPatternFilter.cpp
 * @author Markus Schmidt
 */
#include "module/connectorPatternFilter.h"
#include "db_config.h"

using namespace libMA;

#ifdef WITH_PYTHON

void exportConnectorPatternFilter( py::module& rxPyModuleId )
{
    // export the ConnectorPatternFilter class
    using DBCon = SQLDB<MySQLConDB>;
    exportModule<ConnectorPatternFilter, std::shared_ptr<_SV_DB<DBCon>>>( rxPyModuleId, "ConnectorPatternFilter" );

} // function
#endif