#define USE_DLL_EXPORT
#include "sql_api.h"

/** @brief Single master object for concurrency synchronization. */
#ifdef _MSC_VER
__declspec(dllexport) 
#endif
SQLDBGlobalSync xSQLDBGlobalSync;
