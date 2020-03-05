#pragma once

#include <mysql_con.h>
#include <db_con_pool.h>
#include <sql_api.h>

/// @brief configure the database connection type that is exported to python graph related classes
using DBCon = PooledSQLDBCon<MySQLConDB>;

/// @brief configure the database connection type that is exported to python for non graph related classes
using DBConSingle = SQLDB<MySQLConDB>;