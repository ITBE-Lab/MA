/**
 * @file smith_waterman.h
 * @brief Implements the smith waterman algorithm.
 * @author Markus Schmidt
 */

#ifndef SMITH_WATERMAN
#define SMITH_WATERMAN

#include "module/module.h"
#include "container/alignment.h"

namespace libMABS
{
    class SMW: public Module
    {
        //overload
        std::shared_ptr<Container> execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - NucSeq
         * - Pack
         */
        ContainerVector getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - ContainerVector(Alignment)
         */
        std::shared_ptr<Container> getOutputType() const;

        std::string getName() const
        {
            return "SMW";
        }
    };//class
}//namespace libMABS

void exportSMW();

#endif