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
template <typename TP_CONTAINER> class Lock : public Module<TP_CONTAINER, false, TP_CONTAINER>
{
  public:
    /**
     * @brief create a new Lock.
     * @details
     * Needs information about the container type it shall lock.
     */
    Lock( )
    {} // constructor

    virtual std::shared_ptr<TP_CONTAINER> EXPORTED execute( std::shared_ptr<TP_CONTAINER> pInput )
    {
        // DEBUG( std::cout << "lock" << std::endl; )
        // locking in the container is done automatically by the pledge
        return pInput;
    } // method

}; // class

/**
 * @brief UnLock for enabeling sub-computational graphs
 * @details
 * @see lock
 */
template <typename TP_CONTAINER> class UnLock : public Module<TP_CONTAINER, true, TP_CONTAINER>
{
  public:
    const std::shared_ptr<BasePledge> pLockPledge;

    /**
     * @brief create a new UnLock.
     * @details
     * Takes the Lock it shall unlock as input.
     */
    UnLock( std::shared_ptr<BasePledge> pLockPledge ) : pLockPledge( pLockPledge )
    {} // constructor

    virtual std::shared_ptr<TP_CONTAINER> EXPORTED execute( std::shared_ptr<TP_CONTAINER> pIn )
    {
        // DEBUG( std::cout << "unlock" << std::endl; )
        // unlock the given lock
        pLockPledge->reset( );

        return pIn;
    } // method
}; // class

/**
 * @brief Get a specific tuple element
 * @details
 * the tuple element must contain shared pointers of type TP_TUPLE::value_type
 * the tuple must implement operator[].
 */
template <typename TP_TUPLE, size_t IDX>
class TupleGet : public Module<typename TP_TUPLE::value_type::element_type, true, TP_TUPLE>
{
  public:
    virtual typename TP_TUPLE::value_type EXPORTED execute( std::shared_ptr<TP_TUPLE> pIn )
    {
        return ( *pIn )[ IDX ];
    } // method
}; // class

} // namespace libMA

#ifdef WITH_PYTHON
void exportSplitter( );
#endif

#endif