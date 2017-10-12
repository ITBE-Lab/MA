/** 
 * @file export.h
 * @brief Provides the boost-python main method.
 * @author Markus Schmidt
 * @details 
 * Calls the export functions of the various Modules and Containers.
 */

/**
 * @mainpage Module -ar LineSweep -ing Aligner using Segmentation
 * @section intro_sec Introduction
 *
 * LAuS is a ....
 *
 * @section install_sec Installation
 * 
 * pip install LAuS
 */

/**
 * @defgroup export
 * @brief functions that are used to export Container and Module classes to Python
 * @details
 * When the library is imported in python we need to tell python which classes functions etc.
 * we provide.
 */

#ifndef EXPORT_H
#define EXPORT_H

#include "getAnchors.h"
#include "needlemanWunsch.h"
#include "linesweep.h"
#include "segmentation.h"
#include "bucketing.h"


#endif