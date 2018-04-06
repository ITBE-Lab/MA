/**
 * @file sw_gpu.h
 * @author Markus Schmidt
 * @brief execute a given module on each element of a given libMA::ContainerVector.
 */

#ifndef SW_GPU_H
#define SW_GPU_H

#include "module/module.h"
#include "container/nucSeq.h"

/**
 * @brief Exposes the SweepAllReturnBest @ref Module "module" to boost python.
 * @ingroup export
 */
void exportSW_GPU();

#endif
