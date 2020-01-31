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
#ifndef USE_NEW_DB_API
    exportModule<ConnectorPatternFilter, std::shared_ptr<SV_DB>>( rxPyModuleId, "ConnectorPatternFilter" );
#else
    using DBCon = SQLDB<MySQLConDB>;
	exportModule<ConnectorPatternFilter, std::shared_ptr<_SV_DB<DBCon>>>(rxPyModuleId, "ConnectorPatternFilter");
#endif

} // function
#endif