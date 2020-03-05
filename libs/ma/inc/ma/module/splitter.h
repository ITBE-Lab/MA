/**
 * @file splitter.h
 * @brief exposes exportSplitter
 */
#pragma once

#include "ms/module/module.h"

#ifdef WITH_PYTHON
void exportSplitter( py::module& rxPyModuleId );
#endif