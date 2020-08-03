/* Authors: Arne Kutzner and Markus Schmidt
 * Created: Jan. 2020
 * This file is part of the ITBE-Lab code collection.
 * MIT License
 * @file mysql.h
 * @brief Support of the MySQL database engine.
 */

#pragma once

 // Inform later imports about the SQL backend
#define WITH_MYSQL

#include "mysql_core.h"  // core code for MySQL engine support
#include "mysql_spatial.h" // spatial datatype support for MySQL

