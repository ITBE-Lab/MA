/**
 * @file execOnVector.h
 * @author Markus Schmidt
 * @brief execute a given module on each element of a given libMA::ContainerVector.
 */

#ifndef EXEC_ON_VECTOR_H
#define EXEC_ON_VECTOR_H

#include "module/module.h"
#include "util/threadPool.h"

namespace libMA
{
/**
 * @brief Execute some Module on all elements of a ContainerVector.
 * @details
 * Replaces the first input of the given Module is with a ContainerVector.
 * Then runs the given Module with each element of the vector and all other inputs.
 *
 * Can optionally sort the elements.
 * Sorts in increasing order -> the last element is the one with the highest score.
 * This is necessary since the output is a vector,
 * thus if one wants to extract the largest n elements theu need to be at the end of the vector.
 */
class ExecOnVec : public Module
{
  private:
    std::shared_ptr<Module> pModule;
    bool sort;
    unsigned int nMany; // 0 means return all

  public:
    ExecOnVec( std::shared_ptr<Module> pModule, bool sort = false, unsigned int nMany = 0 )
        : Module( ), pModule( pModule ), sort( sort ), nMany( nMany )
    {} // constructor

    // overload
    std::shared_ptr<Container> EXPORTED execute( std::shared_ptr<ContainerVector> vpInput );

    /**
     * @brief Used to check the input of execute.
     * @details
     * Returns:
     * - ContainerVector(<first input of given Module>)
     * - \<second input of given Module\>
     * - ...
     * - \<last input of given Module\>
     */
    ContainerVector EXPORTED getInputType( ) const;

    /**
     * @brief Used to check the output of execute.
     * @details
     * Returns:
     * - ContainerVector(<output of given Module>)
     */
    std::shared_ptr<Container> EXPORTED getOutputType( ) const;

    std::string getName( ) const
    {
        return "ExecOnVec(" + pModule->getName( ) + ")";
    }

    std::string getFullDesc( ) const
    {
        return std::string( "ExecOnVec(" ) + pModule->getFullDesc( ) + "," +
               std::to_string( sort ) + "," + std::to_string( nMany ) + ")";
    } // function
}; // class

/**
 * @brief Used to extract the last element of a ContainerVector.
 */
class Tail : public Module
{
  private:
    std::shared_ptr<Container> type;

  public:
    Tail( std::shared_ptr<Container> type ) : Module( ), type( type )
    {} // constructor

    // overload
    std::shared_ptr<Container> EXPORTED execute( std::shared_ptr<ContainerVector> vpInput );

    /**
     * @brief Used to check the input of execute.
     * @details
     * Returns:
     * - ContainerVector(x)
     */
    ContainerVector EXPORTED getInputType( ) const;

    /**
     * @brief Used to check the output of execute.
     * @details
     * Returns:
     * - x
     */
    std::shared_ptr<Container> EXPORTED getOutputType( ) const;

    std::string getName( ) const
    {
        return "Tail";
    }
}; // class
} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief Exposes the SweepAllReturnBest @ref Module "module" to boost python.
 * @ingroup export
 */
void exportExecOnVector( );
#endif

#endif
