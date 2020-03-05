/**
 * @file splitter.h
 * @brief exposes exportSplitter
 */
#pragma once

#ifdef WITH_PYTHON
void exportSplitter( py::module& rxPyModuleId );
#endif