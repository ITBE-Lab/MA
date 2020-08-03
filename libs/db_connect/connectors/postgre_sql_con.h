/* Authors: Arne Kutzner and Markus Schmidt
 * Created: July. 2020
 * This file is part of the ITBE-Lab code collection.
 * MIT License
 * @file mysql.h
 * @brief Support of the PostgreSQL database engine.
 */

#pragma once

 // Inform later imports about the SQL backend
#define POSTGRESQL

#include "postgre_sql_core.h"  // core code for PostgreSQL engine support
#include "postgre_sql_spatial.h" // spatial datatype support for PostgreSQL
