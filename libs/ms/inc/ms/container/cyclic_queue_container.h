/**
 * @file cyclic_queue_container.h
 * @brief Implements a queue of containers.
 * @author Markus Schmidt
 */
#pragma once

#include "ms/container/container.h"
#include <condition_variable>
#include <iomanip>
#include <mutex>
#include <queue>

namespace libMS
{

/** @brief a queue that cycles it's elements between push and pop operations.
 *  @details
 *  There are two goals:
 *  - provide the threads with separate queue elements with as little locking as possible.
 *  - don't use more elements than necessary
 * To achieve this, there are two queues in this class:
 * - a priority queue for elements that have already been touched
 * - a low priority queue for elements that are un touched
 */
template <class ContentType> class CyclicQueueBase : public Container
{
    /** @brief low priority queue */
    std::shared_ptr<ContainerVector<std::shared_ptr<ContentType>>> pUntouchedContainers;
    /** @brief used for this->iter() */
    std::shared_ptr<ContainerVector<std::shared_ptr<ContentType>>> pAllContainers;

    std::mutex xMutex;
    std::condition_variable xCondition; // For wait, notify synchronization purposes.

    /** @brief high priority queue */
    std::queue<std::shared_ptr<ContentType>> xTouchedContainers;

  public:
    /** @brief number of elements that have been removed from the queue but not returned */
    size_t uiUnfinishedContainers;

    /** @brief construct a cyclic queue
     * @details
     * Inserts all elements in pContent into the low priority queue
     */
    CyclicQueueBase( std::shared_ptr<ContainerVector<std::shared_ptr<ContentType>>> pContent )
        : pUntouchedContainers( pContent ), pAllContainers( pContent ), uiUnfinishedContainers( pContent->size( ) )
    {} // constructor
    /** @brief construct a cyclic queue
     * @details
     * Both queues are empty in this case
     */
    CyclicQueueBase( )
        : pUntouchedContainers( std::make_shared<ContainerVector<std::shared_ptr<ContentType>>>( ) ),
          pAllContainers( std::make_shared<ContainerVector<std::shared_ptr<ContentType>>>( ) ),
          uiUnfinishedContainers( 0 )
    {} // constructor

    /** @brief removed and returns the head of the queue
     * @details
     * If the priority queue is not empty it's head is removed and returned.
     * Otherwise the head of the low priority queue is removed and returned.
     * If both queues are empty an exception is thrown
     */
    inline std::shared_ptr<ContentType> pop( )
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

    /** @brief puts an element into the queue
     * @details
     * The element is appended to the priority queue
     * @note this can only be called with elements that are taken from this queue via pop()
     */
    inline void push( std::shared_ptr<ContentType> pContainer )
    {
        std::unique_lock<std::mutex> xLock( xMutex );
        xTouchedContainers.push( pContainer );
        xCondition.notify_one( );
    } // method

    /** @brief puts an element into the queue
     * @details
     * The element is appended to the low priority queue.
     * @note this cannot be called with elements that have been taken from this queue via pop()
     */
    inline void add( std::shared_ptr<ContentType> pContainer )
    {
        std::unique_lock<std::mutex> xLock( xMutex );
        pUntouchedContainers->push_back( pContainer );
        pAllContainers->push_back( pContainer );
        uiUnfinishedContainers++;
        xCondition.notify_one( );
    } // method

    inline void informThatContainerIsFinished( )
    {
        std::unique_lock<std::mutex> xLock( xMutex );
        uiUnfinishedContainers--;
        if( uiUnfinishedContainers == 0 )
            xCondition.notify_all( );
    } // method

    inline size_t numTouched( )
    {
        std::unique_lock<std::mutex> xLock( xMutex );
        return xTouchedContainers.size( );
    } // method

    inline size_t numUnTouched( )
    {
        std::unique_lock<std::mutex> xLock( xMutex );
        return pUntouchedContainers->size( );
    } // method

    inline size_t numUnifinished( )
    {
        std::unique_lock<std::mutex> xLock( xMutex );
        return uiUnfinishedContainers;
    } // method

    template <typename F> inline void iterTouched( F&& func )
    {
        std::unique_lock<std::mutex> xLock( xMutex );
        for( size_t uiI = 0; uiI < xTouchedContainers.size( ); uiI++ )
        {
            func( xTouchedContainers.front( ) );
            // there is no regular iteration through the queue so we cycle all elements
            xTouchedContainers.push( xTouchedContainers.front( ) );
            xTouchedContainers.pop( );
        } // for
    } // method
    template <typename F> inline void iter( F&& func )
    {
        std::unique_lock<std::mutex> xLock( xMutex );
        for( auto pContainer : *pAllContainers )
            func( pContainer );
    } // method

}; // class

// switch between profiled and un-profiled version of cyclic queue
#if 0 // 1-> profiled / 0 -> un-profiled

/** @brief see CyclicQueueBase for documentation
 *  @details
 *  The profiled CyclicQueue. This class implements a profiler for the queue.
 */
template <class ContentType> class CyclicQueue : public Container
{
    using duration = std::chrono::duration<double>;
    using time_point = std::chrono::time_point<std::chrono::steady_clock, duration>;

    CyclicQueueBase<ContentType> xQueue;

    size_t uiActiveContainers = 0;
    size_t numTouchedLast = 0;
    size_t numUnTouchedLast = 0;
    std::mutex xMutex;
    time_point xLastTime;
    time_point xInitTimePoint;
    std::vector<duration> vTimeActive;
    std::vector<duration> vTimeTouched;
    std::vector<duration> vTimeUnTouched;
    duration xPopTimeTotal;
    duration xPushTimeTotal;

    /**
     * @brief updates the time counters
     * @details
     * Has to be called immediately after a change to xQueue, but before uiActiveContainers is adjusted.
     * This is necessary, since xQueue.pop() triggers a wait operation within the queue,
     * so calling updateTimes just before xQueue operations is not possible.
     */
    inline void updateTimes( )
    {
        auto xNow = std::chrono::steady_clock::now( );

        while( vTimeActive.size( ) <= uiActiveContainers )
            vTimeActive.push_back( duration::zero( ) );
        vTimeActive[ uiActiveContainers ] += xNow - xLastTime;

        auto uiI = numTouchedLast + uiActiveContainers;
        while( vTimeTouched.size( ) <= uiI )
            vTimeTouched.push_back( duration::zero( ) );
        vTimeTouched[ uiI ] += xNow - xLastTime;

        uiI = numUnTouchedLast;
        while( vTimeUnTouched.size( ) <= uiI )
            vTimeUnTouched.push_back( duration::zero( ) );
        vTimeUnTouched[ uiI ] += xNow - xLastTime;

        xLastTime = xNow;
        numTouchedLast = xQueue.numTouched( );
        numUnTouchedLast = xQueue.numUnTouched( );
    } // method

    inline void initTime( )
    {
        xLastTime = std::chrono::steady_clock::now( );
        xInitTimePoint = xLastTime;
        vTimeActive.clear( );
        vTimeTouched.clear( );
        vTimeUnTouched.clear( );

        xPopTimeTotal = duration::zero( );
        xPushTimeTotal = duration::zero( );
    } // method

    inline void printArray( std::ostream& xOut, const std::string& sName, const std::vector<duration>& vArr )
    {
        duration xTotal = duration::zero( );
        for( size_t uiI = 0; uiI < vArr.size( ); uiI++ )
            xTotal += vArr[ uiI ];
        xOut << sName << "\t";
        for( size_t uiI = 0; uiI < vArr.size( ); uiI++ )
        {
            double dPercentage = 100.0 * vArr[ uiI ].count( ) / xTotal.count( );
            if( dPercentage > 5.0 )
                xOut << uiI << ": " << dPercentage << "%\t";
        } // for
        xOut << std::endl;
    } // method

    inline void printTime( std::ostream& xOut )
    {
        xOut << "CyclicQueue profiler:" << std::endl;
        xOut << std::fixed << std::setprecision( 2 );
        std::cout << "Time spent in (accumulated over threads): push(): "
                  << 100.0 * xPushTimeTotal.count( ) / ( xLastTime - xInitTimePoint ).count( )
                  << "% pop(): " << 100.0 * xPopTimeTotal.count( ) / ( xLastTime - xInitTimePoint ).count( ) << "%"
                  << std::endl;
        // two digits precision double output
        printArray( xOut, "#Active   ", vTimeActive );
        printArray( xOut, "#Touched  ", vTimeTouched );
        printArray( xOut, "#UnTouched", vTimeUnTouched );
        xOut << std::defaultfloat;
    } // method

    inline void top( std::ostream& xOut )
    {
        if( ( xLastTime - xInitTimePoint ).count( ) > 10 /*amount of seconds between prints*/ )
        {
            printTime( xOut );
            initTime( );
        } // if
    } // method

  public:
    CyclicQueue( std::shared_ptr<ContainerVector<std::shared_ptr<ContentType>>> pContent ) : xQueue( pContent )
    {
        numTouchedLast = xQueue.numTouched( );
        numUnTouchedLast = xQueue.numUnTouched( );
        initTime( );
    } // constructor

    CyclicQueue( ) : xQueue( )
    {
        numTouchedLast = xQueue.numTouched( );
        numUnTouchedLast = xQueue.numUnTouched( );
        initTime( );
    } // constructor

    ~CyclicQueue( )
    {
        printTime( std::cout );
    } // destructor

    inline std::shared_ptr<ContentType> pop( )
    {
        auto xStart = std::chrono::steady_clock::now( );
        auto pRet = xQueue.pop( );
        auto xDur = std::chrono::steady_clock::now( ) - xStart;
        if( pRet != nullptr )
        {
            std::unique_lock<std::mutex> xLock( xMutex );
            xPopTimeTotal += xDur;
            updateTimes( );
            uiActiveContainers++;
        } // if
        return pRet;
    } // method

    inline void push( std::shared_ptr<ContentType> pContainer )
    {
        auto xStart = std::chrono::steady_clock::now( );
        xQueue.push( pContainer );
        auto xDur = std::chrono::steady_clock::now( ) - xStart;
        std::unique_lock<std::mutex> xLock( xMutex );
        xPushTimeTotal += xDur;
        updateTimes( );
        top( std::cout );
        assert( uiActiveContainers > 0 );
        uiActiveContainers--;
    } // method

    inline void add( std::shared_ptr<ContentType> pContainer )
    {
        xQueue.add( pContainer );
    } // method

    inline void informThatContainerIsFinished( )
    {
        auto xStart = std::chrono::steady_clock::now( );
        xQueue.informThatContainerIsFinished( );
        auto xDur = std::chrono::steady_clock::now( ) - xStart;
        std::unique_lock<std::mutex> xLock( xMutex );
        xPushTimeTotal += xDur;
        updateTimes( );
        assert( uiActiveContainers > 0 );
        uiActiveContainers--;
    } // method

    inline size_t numTouched( )
    {
        return xQueue.numTouched( );
    } // method

    inline size_t numUnTouched( )
    {
        return xQueue.numUnTouched( );
    } // method

    inline size_t numUnifinished( )
    {
        return xQueue.numUnifinished( );
    } // method
    template<typename F>
    inline void iterTouched( F&& func )
    {
        xQueue.iterTouched( func );
    } // method
    template<typename F>
    inline void iter( F&& func )
    {
        xQueue.iter( func );
    } // method
}; // class

#else

/** @brief see CyclicQueueBase for documentation
 *  @details
 *  The un-profiled CyclicQueue. In this cases no changes to CyclicQueueBase are needed.
 */
template <class ContentType> using CyclicQueue = CyclicQueueBase<ContentType>;

#endif

} // namespace libMS