/**
 * @file sw_gpu.h
 * @author Markus Schmidt
 * @brief execute a given module on each element of a given libMA::ContainerVector.
 */

#ifndef SW_GPU_H
#define SW_GPU_H

#include "module/module.h"
#include "container/nucSeq.h"

namespace libMA
{
    class SW_GPU: public Module
    {
    public:
        SW_GPU()
                :
            Module()
        {}//constructor

        //overload
        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> vpInput);
        
        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - ContainerVector(x)
         */
        ContainerVector EXPORTED getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - x
         */
        std::shared_ptr<Container> EXPORTED getOutputType() const;
        
        std::string getName() const
        {
            return "SW_GPU";
        }
    };//class
}//namespace libMA

/**
 * @brief Exposes the SweepAllReturnBest @ref Module "module" to boost python.
 * @ingroup export
 */
void exportSW_GPU();

#endif
