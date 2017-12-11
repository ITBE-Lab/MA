#ifndef SMITH_WATERMAN
#define SMITH_WATERMAN

#include "cppModule.h"
#include "alignment.h"

void exportSMW();

class SMW: public CppModule
{
    //overload
    std::shared_ptr<Container> execute(ContainerVector vpInput);

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

#endif