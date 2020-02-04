#pragma once

#include <db_con_pool.h>
#include <common.h>

/// @brief configure the database connection type that is exported to python graph related classes
using DBCon = PooledSQLDBCon<MySQLConDB>;
/// @brief configure the database connection type that is exported to python for non graph related classes
using DBConSingle = SQLDB<MySQLConDB>;