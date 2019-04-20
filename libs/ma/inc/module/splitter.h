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
    Lock( const ParameterSetManager& rParameters )
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
    UnLock( const ParameterSetManager& rParameters, std::shared_ptr<BasePledge> pLockPledge )
        : pLockPledge( pLockPledge )
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
class TupleGet : public Module<typename TP_TUPLE::value_type::element_type, false, TP_TUPLE>
{
  public:
    TupleGet( const ParameterSetManager& rParameters )
    {} // constructor

    virtual typename TP_TUPLE::value_type EXPORTED execute( std::shared_ptr<TP_TUPLE> pIn )
    {
        return ( *pIn )[ IDX ];
    } // method
}; // class

/**
 * @brief Split a ContainerVector into its elements
 * @details
 */
template <typename TP> class Splitter : public Module<TP, true, ContainerVector<std::shared_ptr<TP>>>
{
  public:
    Splitter( const ParameterSetManager& rParameters )
    {} // constructor

    // @override
    virtual bool requiresLock( ) const
    {
        return true;
    } // function

    virtual typename std::shared_ptr<TP> EXPORTED execute( std::shared_ptr<ContainerVector<std::shared_ptr<TP>>> pIn )
    {
        if( pIn->empty( ) )
            // if we reach this point we have read all content vector
            throw AnnotatedException( "Tried to extract element from empty vector" );
        auto pBack = pIn->back( );
        pIn->pop_back( );
        if( pIn->empty( ) )
            this->setFinished( );
        return pBack;
    } // method
}; // class

/**
 * @brief Split a ContainerVector into its elements
 * @details
 */
template <typename TP> class StaticSplitter : public Module<TP, true>
{
  public:
    std::shared_ptr<ContainerVector<std::shared_ptr<TP>>> pIn;
    StaticSplitter( const ParameterSetManager& rParameters, std::shared_ptr<ContainerVector<std::shared_ptr<TP>>> pIn )
        : pIn( pIn )
    {} // constructor

    // @override
    virtual bool requiresLock( ) const
    {
        return true;
    } // function

    typename std::shared_ptr<TP> execute( void )
    {
        if( pIn->empty( ) )
            // if we reach this point we have read all content vector
            throw AnnotatedException( "Tried to extract element from empty vector" );
        auto pBack = pIn->back( );
        pIn->pop_back( );
        if( pIn->empty( ) )
            this->setFinished( );
        return pBack;
    } // method
}; // class

/**
 * @brief Get a specific tuple element
 * @details
 * the tuple element must contain shared pointers of type TP_TUPLE::value_type
 * the tuple must implement operator[].
 */
template <typename... TP_VEC_CONTENT> class Collector : public Module<Container, false, TP_VEC_CONTENT...>
{
  public:
    std::vector<std::tuple<std::shared_ptr<TP_VEC_CONTENT>...>> vCollection;
    std::shared_ptr<std::mutex> pMutex;

    Collector( const ParameterSetManager& rParameters ) : pMutex( new std::mutex )
    {} // constructor

    virtual std::shared_ptr<Container> execute( std::shared_ptr<TP_VEC_CONTENT>... pIn )
    {
        std::lock_guard<std::mutex> xGuard( *pMutex );
        vCollection.push_back( std::make_tuple( pIn... ) );
        return std::make_shared<Container>( );
    } // method
}; // class

/**
 * @brief Joins two ends of computational graphs
 * @details
 * always returns an empty container
 */
template <typename... TP_VEC_CONTENT> class Join : public Module<Container, false, TP_VEC_CONTENT...>
{
  public:
    Join( const ParameterSetManager& rParameters )
    {} // constructor

    virtual typename std::shared_ptr<Container> EXPORTED execute( std::shared_ptr<TP_VEC_CONTENT>... pIn )
    {
        return std::make_shared<Container>( );
    } // method
}; // class


} // namespace libMA

#ifdef WITH_PYTHON
#ifdef WITH_BOOST
void exportSplitter( );
#else
void exportSplitter( py::module& rxPyModuleId );
#endif
#endif

#endif