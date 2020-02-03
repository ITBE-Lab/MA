#pragma once

#include <db_con_pool.h>

/// @brief configure the database connection type that is exported to python
using DBCon = PooledSQLDBCon<MySQLConDB>;