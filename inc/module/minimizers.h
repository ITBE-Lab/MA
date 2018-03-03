/** 
 * @file binarySeeding.h
 * @brief Implements a segmentation algorithm.
 * @author Markus Schmidt
 */
#ifndef MINIMIZERS_H
#define MINIMIZERS_H

#include "util/system.h"
#include "module/module.h"
#include "container/segment.h"
#include "container/minimizersHash.h"
#include "util/threadPool.h"

namespace libMA
{

    class Minimizers : public Module{
    public:
        const static unsigned int w = 5, k = 15;
        bool bPrint = false;

        Minimizers()
        {}//constructor
        
        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - NucSeq
         */
        ContainerVector EXPORTED getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - SegmentVector
         */
        std::shared_ptr<Container> EXPORTED getOutputType() const;

        std::string getName() const
        {
            return "Minimizers";
        }
    };//class

    class MinimizersToSeeds : public Module{
    public:

        MinimizersToSeeds()
        {}//constructor
        
        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - FMIndex
         * - NucSeq
         */
        ContainerVector EXPORTED getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - SegmentVector
         */
        std::shared_ptr<Container> EXPORTED getOutputType() const;

        std::string getName() const
        {
            return "MinimizersToSeeds";
        }
    };//class

}//namespace


/**
 * @brief exports the Segmentation @ref Module "module" to python.
 * @ingroup export
 */
void exportMinimizers();

#endif