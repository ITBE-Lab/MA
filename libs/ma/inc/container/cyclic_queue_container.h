/**
 * @file cyclic_queue_container.h
 * @brief Implements a queue of containers.
 * @author Markus Schmidt
 */
#pragma once

#include "container/container.h"
#include <mutex>
#include <queue>
#include <condition_variable>

namespace libMA
{

template <class ContentType> class CyclicQueue : public Container
{
    std::shared_ptr<ContainerVector<std::shared_ptr<ContentType>>> pUntouchedContainers;

    std::mutex xMutex;
    std::condition_variable xCondition; // For wait, notify synchronization purposes.
    std::queue<std::shared_ptr<ContentType>> xTouchedContainers;

  public:
    size_t uiUnfinishedContainers;

    CyclicQueue( std::shared_ptr<ContainerVector<std::shared_ptr<ContentType>>> pContent )
        : pUntouchedContainers( pContent ), uiUnfinishedContainers( pContent->size( ) )
    {} // constructor

    CyclicQueue( )
        : pUntouchedContainers( std::make_shared<ContainerVector<std::shared_ptr<ContentType>>>( ) ),
          uiUnfinishedContainers( 0 )
    {} // constructor

    std::shared_ptr<ContentType> pop( )
    {
        std::unique_lock<std::mutex> xLock( xMutex );

        while( xTouchedContainers.empty( ) && pUntouchedContainers->empty( ) && uiUnfinishedContainers > 0 )
            xCondition.wait( xLock );
        if( uiUnfinishedContainers == 0 ) // we are completely done
            return nullptr;

        if( !xTouchedContainers.empty( ) )
        {
            auto pRet = xTouchedContainers.front( );
            xTouchedContainers.pop( );
            return pRet;
        } // if
        else if( !pUntouchedContainers->empty( ) )
        {
            auto pRet = pUntouchedContainers->back( );
            pUntouchedContainers->pop_back( );
            return pRet;
        } // else if

        throw std::runtime_error( "should never be able to read this line" );
    } // method

    void push( std::shared_ptr<ContentType> pContainer )
    {
        std::unique_lock<std::mutex> xLock( xMutex );
        xTouchedContainers.push( pContainer );
        xCondition.notify_one( );
    } // method

    void add( std::shared_ptr<ContentType> pContainer )
    {
        std::unique_lock<std::mutex> xLock( xMutex );
        pUntouchedContainers->push_back( pContainer );
        uiUnfinishedContainers++;
        xCondition.notify_one( );
    } // method

    void informThatContainerIsFinished( )
    {
        std::unique_lock<std::mutex> xLock( xMutex );
        uiUnfinishedContainers--;
        if( uiUnfinishedContainers == 0 )
            xCondition.notify_all( );
    } // method

}; // class

} // namespace libMA