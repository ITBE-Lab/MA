/**
 * @file svDb.h
 * @details
 * The database interface for the structural variant caller
 */
#pragma once

#include "module/module.h"

#ifdef WITH_PYTHON
void exportSoCDbWriter( py::module& rxPyModuleId );
#endif
