/* Authors: Arne Kutzner and Markus Schmidt
 * Created: July 2020
 * This file is part of the ITBE-Lab code collection.
 * MIT License
 * @db_base.h
 * @brief This file must always be included first. It selects the database backend.
 */
#pragma once

//#define USE_MYSQL
#define USE_PG

#ifdef USE_PG
#include "postgre_sql_con.h"
using DBConn = PostgreSQLDBCon;
#endif

#ifdef USE_MYSQL
#include <mysql_con.h>
using DBConn = MySQLConDB;
#endif

#define DB_BASE_INCLUDED