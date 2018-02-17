/**
 * @file smith_waterman.h
 * @brief Implements the smith waterman algorithm.
 * @author Markus Schmidt
 */

#ifndef SMITH_WATERMAN
#define SMITH_WATERMAN

#include "module/module.h"
#include "container/alignment.h"

namespace libMA
{
    class SMW: public Module
    {
    public:
        bool bBacktrack = true;
        DEBUG(
            bool bPrint = false;
        )//DEBUG

        SMW(bool bBacktrack)
                :
            bBacktrack(bBacktrack)
        {}//constructor

        //overload
        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - NucSeq
         * - Pack
         */
        ContainerVector EXPORTED getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - ContainerVector(Alignment)
         */
        std::shared_ptr<Container> EXPORTED getOutputType() const;

        std::string getName() const
        {
            return "SMW";
        }
    };//class
}//namespace libMA

void exportSMW();

#endif