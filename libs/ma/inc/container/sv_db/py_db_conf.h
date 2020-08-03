#pragma once

#include <db_base.h>
#include <db_con_pool.h>
// #include <sql_api.h>

/// @brief configure the database connection type that is exported to python graph related classes
using DBCon = PooledSQLDBCon<DBConn>;
/// @brief configure the database connection type that is exported to python for non graph related classes
using DBConSingle = SQLDB<DBConn>;