/** 
 * @file needlemanWunsch.h
 * @brief Implements NMW
 * @author Markus Schmidt
 */
#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

#include "alignment.h"
#include "cppModule.h"

/**
 * @brief implements NMW
 * @details
 * Returns a finished alignment if given a sound selection of seeds.
 * @ingroup module
 */
class NeedlemanWunsch : public CppModule
{
public:

    NeedlemanWunsch()
    {}//constructor

    //overload
    std::shared_ptr<Container> execute(ContainerVector vpInput);

    /**
     * @brief Used to check the input of execute.
     * @details
     * Returns:
     * - Seeds
     * - NucSeq
     * - Pack
     */
    ContainerVector getInputType() const;

	/**
	 * @brief Used to check the output of execute.
	 * @details
	 * Returns:
	 * - Alignment
	 */
    std::shared_ptr<Container> getOutputType() const;

    std::string getName() const
    {
        return "NeedlemanWunsch";
    }
};//class

/**
 * @brief Exposes the Alignment container to boost python.
 * @ingroup export
 */
void exportNeedlemanWunsch();

#endif