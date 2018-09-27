/**
 * @file splitter.h
 * @brief implements the ability to make sub-sections in computational graphs
 * @author Markus Schmidt
 * @todo there should be an actual sub comp. graph class instead of
 * exposing locks, unlocks splitters and so on...
 */
#ifndef SPLITTER_H
#define SPLITTER_H

#include "container/nucSeq.h"
#include "module/module.h"
#include <fstream>
#include <mutex>

namespace libMA
{
/**
 * @brief splits the a container vector into its elements
 * @details
 * Intended usage:
 * Given a container vectors of queries, this returns one query each time it is executed.
 * This way multiple sub comp. graphs can take queries from the same vector and align them
 * simultaneously.
 */
class Splitter : public Module
{
  public:
    std::shared_ptr<Pledge> pVec;
    std::shared_ptr<std::mutex> pMutex;

    /**
     * @brief create a new splitter.
     * @details
     * Takes the vector that shall be split as input.
     */
    Splitter( std::shared_ptr<Pledge> pVec ) : pVec( pVec ), pMutex( new std::mutex )
    {} // constructor

    std::shared_ptr<Container> EXPORTED execute( std::shared_ptr<ContainerVector> vpInput );

    /**
     * @brief Used to check the input of execute.
     * @details
     * Returns:
     * - Nil
     * The input is given in the contructor...
     */
    ContainerVector EXPORTED getInputType( ) const;

    /**
     * @brief Used to check the output of execute.
     * @details
     * Returns:
     * - ContainerVector(NucSeq)
     */
    std::shared_ptr<Container> EXPORTED getOutputType( ) const;

    std::string getName( ) const
    {
        return "Splitter";
    } // function

    bool outputsVolatile( ) const
    {
        return true;
    } // function

    std::string getFullDesc( ) const
    {
        return std::string( "Splitter" );
    } // function

}; // class


class ReadSplitter : public Module
{
  public:
    ReadSplitter( )
    {} // constructor

    std::shared_ptr<Container> EXPORTED execute( std::shared_ptr<ContainerVector> vpInput );

    /**
     * @brief Used to check the input of execute.
     * @details
     * Returns:
     * - Nil
     * The input is given in the contructor...
     */
    ContainerVector EXPORTED getInputType( ) const;

    /**
     * @brief Used to check the output of execute.
     * @details
     * Returns:
     * - ContainerVector(NucSeq)
     */
    std::shared_ptr<Container> EXPORTED getOutputType( ) const;

    std::string getName( ) const
    {
        return "ReadSplitter";
    } // function

    bool outputsVolatile( ) const
    {
        return true;
    } // function

    std::string getFullDesc( ) const
    {
        return std::string( "ReadSplitter" );
    } // function

}; // class

/**
 * @brief collects the outputs of multiple sub-computational graphs and combines them
 * @details
 * Intended usage:
 * If multiple computational graphs are used to align simultaneously this module collects
 * the output in one vector.
 */
class Collector : public Module
{
  public:
    std::shared_ptr<ContainerVector> pVec;
    std::shared_ptr<std::mutex> pLock;

    /**
     * @brief create a new Collector.
     * @details
     * Needs information about the container type it shall collect.
     */
    Collector( std::shared_ptr<Container> pType )
        : pVec( new ContainerVector( pType ) ), pLock( new std::mutex )
    {} // constructor

    std::shared_ptr<Container> EXPORTED execute( std::shared_ptr<ContainerVector> vpInput );

    /**
     * @brief Used to check the input of execute.
     * @details
     * Returns:
     * - Nil
     * The output is in pVec.
     */
    ContainerVector EXPORTED getInputType( ) const;

    /**
     * @brief Used to check the output of execute.
     * @details
     * Returns:
     * - ContainerVector(NucSeq)
     */
    std::shared_ptr<Container> EXPORTED getOutputType( ) const;

    std::string getName( ) const
    {
        return "Collector";
    } // function

    std::string getFullDesc( ) const
    {
        return std::string( "Collector" );
    } // function

}; // class

/**
 * @brief Lock for enabeling sub-computational graphs
 * @details
 * This gets input from another module and locks this input until the respective UnLock tells
 * it to replace the output with another.
 * Always use in combination with UnLock!
 * Intended usage:
 * A Container element need to be used multiple time in some comp-subgraph.
 * However, each time get() is called to obtain the element a new element is given since
 * the element comes from a volatile Module.
 * This locks this element and therefore enables is re usage.
 * Once the element does not need to be used anymore a UnLock can be used to unlock the element.
 */
class Lock : public Module
{
  public:
    std::shared_ptr<Container> pType;

    /**
     * @brief create a new Lock.
     * @details
     * Needs information about the container type it shall lock.
     */
    Lock( std::shared_ptr<Container> pType ) : pType( pType )
    {} // constructor

    std::shared_ptr<Container> EXPORTED execute( std::shared_ptr<ContainerVector> vpInput );

    /**
     * @brief Used to check the input of execute.
     * @details
     * Returns:
     * - Nil
     */
    ContainerVector EXPORTED getInputType( ) const;

    /**
     * @brief Used to check the output of execute.
     * @details
     * Returns:
     * - ContainerVector(NucSeq)
     */
    std::shared_ptr<Container> EXPORTED getOutputType( ) const;

    std::string getName( ) const
    {
        return "Lock";
    } // function

    std::string getFullDesc( ) const
    {
        return std::string( "Lock" );
    } // function

    // @override
    bool requiresLock( ) const
    {
        return true;
    } // function

}; // class

/**
 * @brief UnLock for enabeling sub-computational graphs
 * @details
 * @see lock
 */
class UnLock : public Module
{
  public:
    std::shared_ptr<Pledge> pLockPledge;

    /**
     * @brief create a new UnLock.
     * @details
     * Takes the Lock it shall unlock as input.
     */
    UnLock( std::shared_ptr<Pledge> pLockPledge ) : pLockPledge( pLockPledge )
    {} // constructor

    std::shared_ptr<Container> EXPORTED execute( std::shared_ptr<ContainerVector> vpInput );

    /**
     * @brief Used to check the input of execute.
     * @details
     * Returns:
     * - Nil
     */
    ContainerVector EXPORTED getInputType( ) const;

    /**
     * @brief Used to check the output of execute.
     * @details
     * Returns:
     * - ContainerVector(NucSeq)
     */
    std::shared_ptr<Container> EXPORTED getOutputType( ) const;

    std::string getName( ) const
    {
        return "UnLock";
    } // function

    std::string getFullDesc( ) const
    {
        return "UnLock";
    } // function

    bool outputsVolatile( ) const
    {
        return true;
    } // function
}; // class

} // namespace libMA

#ifdef WITH_PYTHON
void exportSplitter( );
#endif

#endif