#include "connectors/postgre_sql_core.h"


/** @brief Single master object for concurrency synchronization. */
PG_GLobalEnv xPG_GLobalEnv;