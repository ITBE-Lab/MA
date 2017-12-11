/**
 * TODO:
 * @file execOnVector.h
 */

#ifndef EXEC_ON_VECTOR_H
#define EXEC_ON_VECTOR_H

#include "cppModule.h"
#include "threadPool.h"

/**
 * @brief Execute some CppModule on all elements of a ContainerVector.
 * @details
 * Replaces the first input of the given CppModule is with a ContainerVector.
 * Then runs the given CppModule with each element of the vector and all other inputs.
 * 
 * Can optionally sort the elements.
 * Sorts in increasing order -> the last element is the one with the highest score.
 * This is necessary since the output is a vector, 
 * thus if one wants to extract the largest n elements theu need to be at the end of the vector.
 */
class ExecOnVec: public CppModule
{
private:
    std::shared_ptr<CppModule> pModule;
    bool sort;
    unsigned int nMany;//0 means return all

public:
    ExecOnVec(std::shared_ptr<CppModule> pModule, bool sort, unsigned int nMany)
            :
        CppModule(),
        pModule(pModule),
        sort(sort),
        nMany(nMany)
    {}//constructor

    //overload
    std::shared_ptr<Container> execute(ContainerVector vpInput);

    /**
     * @brief Used to check the input of execute.
     * @details
     * Returns:
     * - ContainerVector(<first input of given Module>)
     * - \<second input of given Module\>
     * - ...
     * - \<last input of given Module\>
     */
    ContainerVector getInputType() const;

	/**
	 * @brief Used to check the output of execute.
	 * @details
	 * Returns:
	 * - ContainerVector(<output of given Module>)
	 */
    std::shared_ptr<Container> getOutputType() const;
    
    std::string getName() const
    {
        return "ExecOnVec( " + pModule->getName() + ")";
    }
};//class

/**
 * @brief Used to extract the last element of a ContainerVector.
 */
class Tail: public CppModule
{
private:
    std::shared_ptr<Container> type;
public:
    Tail(std::shared_ptr<Container> type)
            :
        CppModule(),
        type(type)
    {}//constructor

    //overload
    std::shared_ptr<Container> execute(ContainerVector vpInput);
    
    /**
     * @brief Used to check the input of execute.
     * @details
     * Returns:
     * - ContainerVector(x)
     */
    ContainerVector getInputType() const;

	/**
	 * @brief Used to check the output of execute.
	 * @details
	 * Returns:
	 * - x
	 */
    std::shared_ptr<Container> getOutputType() const;
    
    std::string getName() const
    {
        return "Tail";
    }
};//class

/**
 * @brief Exposes the SweepAllReturnBest @ref CppModule "module" to boost python.
 * @ingroup export
 */
void exportExecOnVector();

#endif
