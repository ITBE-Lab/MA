/**
 * @file connectorPatternFilter.cpp
 * @author Markus Schmidt
 */
#include "msv/module/connectorPatternFilter.h"

using namespace libMSV;
using namespace libMS;

#ifdef WITH_PYTHON

#include "ms/container/sv_db/py_db_conf.h"

void exportConnectorPatternFilter( libMS::SubmoduleOrganizer& xOrganizer )
{
    // export the ConnectorPatternFilter class
    exportModule<ConnectorPatternFilter<DBCon>>( xOrganizer, "ConnectorPatternFilter" );

} // function
#endif