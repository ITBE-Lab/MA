#define USE_DLL_EXPORT
#include "connectors/postgre_sql_core.h"
#include "connectors/postgre_sql_spatial.h"

/** @brief Single master object for concurrency synchronization. */
#ifdef _MSC_VER
__declspec(dllexport) PG_GLobalEnv xPG_GLobalEnv;
#else
PG_GLobalEnv xPG_GLobalEnv;
#endif
