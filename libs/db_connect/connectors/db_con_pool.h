/* Authors: Arne Kutzner and Markus Schmidt
 * Created: Jan. 2020
 * This file is part of the ITBE-Lab code collection.
 * MIT License
 * @file db_con_pool.h
 * @brief Connection pool for databases
 */

#pragma once

#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

#include "sql_api.h" // general API for SQL like DB


// Rename: PooledSQLDBCon => PooledSQLDBCon
template <typename DBImpl> class PooledSQLDBCon : public SQLDB<DBImpl>
{
  private:
    std::shared_ptr<std::mutex> pPoolLock; // the pool lock is shared with the pool and all its connections

  public:
    using DBImplForwarded = DBImpl;

    size_t uiThreadId; // id of the task belonging to the PooledSQLDBCon
    PooledSQLDBCon( const PooledSQLDBCon<DBImpl>& ) = delete; // no copies of pooled connections

    PooledSQLDBCon( std::shared_ptr<std::mutex> pPoolLock, // pass shared pointer to pool mutex
                    const json& jDBConData = json{} ) // FIXME: Pass JSON here
        : SQLDB<DBImpl>( jDBConData ), pPoolLock( pPoolLock )
    {} // constructor

    /** @brief The function func is executed protected by a global lock. */
    template <typename F> void doPoolSafe( F&& func )
    {
        std::lock_guard<std::mutex> lock( *pPoolLock );
        // DEBUG: std::cout << "doPoolSafe in PooledSQLDBCon ..." << std::endl;
        func( );
    } // method

    inline size_t getTaskId( )
    {
        return this->uiThreadId;
    } // method
}; // class


/* Basic characteristics:
 * - The pool synchronizes threads that ask for access, but the pool itself does not comprises threads.
 *   (Therefore, e.g. an additional thread pool is required for delivering multiple thread)
 * Important notice: If a pooled task throws an exception, this exception belongs to the future and must be caught via
 * the future. For more details see:
 * https://stackoverflow.com/questions/16344663/is-there-a-packaged-taskset-exception-equivalent
 *
 */
template <typename DBImpl> class SQLDBConPool
{
#define DO_PROFILE 0
  private:
    using TaskType = std::function<void( std::shared_ptr<PooledSQLDBCon<DBImpl>> )>; // type of the tasks to be executed


    std::vector<std::thread> vWorkers; // worker threads; each thread manages one connection
    std::vector<std::shared_ptr<PooledSQLDBCon<DBImpl>>> vConPool; // pointers to actual connections
    std::vector<std::queue<TaskType>> vqTasks; // queue of lambda functions that shall to executed
#if DO_PROFILE == 1
    using duration = std::chrono::duration<double>;
    using time_point = std::chrono::time_point<std::chrono::steady_clock, duration>;

    std::vector<std::array<duration, 4>> vTimeAmountTasksInQueue;
    std::vector<std::tuple<time_point, duration, bool>> vTimeThreadSleeping;
    std::vector<uint64_t> vAmountTasks;
    std::vector<uint64_t> vAmountBookedTasks;

    time_point xLastTimePoint;
    time_point xInitTimePoint;

    /**
     * @brief initializes the usage measurement
     * @details
     * call this before you start using any of the usage measurement functions
     * calling this again will reset the usage measurement
     *
     * all the Time functions are not threadsave and assume that they are called while the pPoolLock is held
     */
    inline void initTime( size_t uiPoolSize )
    {
        vTimeAmountTasksInQueue.clear( );
        vAmountTasks.clear( );
        xLastTimePoint = std::chrono::steady_clock::now( );
        xInitTimePoint = xLastTimePoint;
        for( size_t uiI = 0; uiI < uiPoolSize; uiI++ )
        {
            vTimeAmountTasksInQueue.emplace_back( );
            // if init is called while threads are sleeping this will
            if( vTimeThreadSleeping.size( ) <= uiI )
                vTimeThreadSleeping.emplace_back( xLastTimePoint, duration::zero( ), false );
            else
            {
                std::get<0>( vTimeThreadSleeping[ uiI ] ) = xLastTimePoint;
                std::get<1>( vTimeThreadSleeping[ uiI ] ) = duration::zero( );
            } // else
            vAmountTasks.emplace_back( 0 );
            if( vAmountBookedTasks.size( ) <= uiI )
                vAmountBookedTasks.emplace_back( 0 );
        } // for
        vAmountTasks.emplace_back( );
        vTimeAmountTasksInQueue.emplace_back( );
    } // method

    /**
     * @brief call this whenever a thread gets a new task
     */
    inline void incTasks( size_t uiThreadId )
    {
        vAmountTasks[ uiThreadId ]++;
    } // method

    /**
     * @brief call this whenever a thread id is booked
     */
    inline void bookTask( size_t uiThreadId )
    {
        vAmountBookedTasks[ uiThreadId ]++;
    } // method

    /**
     * @brief call this just before the the size of any task queue changes
     */
    inline void updateTime( )
    {
        auto xNow = std::chrono::steady_clock::now( );
        auto xTimePassed = xNow - xLastTimePoint;

        for( size_t uiI = 0; uiI < vTimeAmountTasksInQueue.size( ); uiI++ )
            vTimeAmountTasksInQueue[ uiI ][ std::min( vqTasks[ uiI ].size( ), (size_t)3 ) ] += xTimePassed;

        xLastTimePoint = xNow;
    } // method

    /**
     * @brief call this just before a thread goes to sleep
     */
    inline void goingToSleep( size_t uiThreadId )
    {
        std::get<0>( vTimeThreadSleeping[ uiThreadId ] ) = std::chrono::steady_clock::now( );
        std::get<2>( vTimeThreadSleeping[ uiThreadId ] ) = true;
    } // method

    /**
     * @brief call this just after a thread wakes up
     */
    inline void wakingUp( size_t uiThreadId )
    {
        auto xNow = std::chrono::steady_clock::now( );
        // sanity check
        std::get<1>( vTimeThreadSleeping[ uiThreadId ] ) += xNow - std::get<0>( vTimeThreadSleeping[ uiThreadId ] );
        std::get<2>( vTimeThreadSleeping[ uiThreadId ] ) = false;
    } // method

    /**
     * @brief call this to print a summary of the activity within the pool
     */
    inline void printTime( std::ostream& xOut )
    {
        xOut << "SQLDBConPool profiler: " << std::endl;
        xOut << "tID\t#tsk/s\t#booked\t%active\t%qlen=0\t%qlen=1\t%qlen=2\t%qlen>2" << std::endl;

        // two digits precision double output
        xOut << std::fixed << std::setprecision( 2 );

        double fTotalTime = ( xLastTimePoint - xInitTimePoint ).count( );
        auto xNow = std::chrono::steady_clock::now( );

        for( size_t uiI = 0; uiI < vTimeAmountTasksInQueue.size( ); uiI++ )
        {
            // print the values for tID #tasks #book %act
            auto fTasksPerSecond = vAmountTasks[ uiI ] / fTotalTime;
            if( uiI < vTimeThreadSleeping.size( ) )
            {
                // get the amount of time thread uiI was active
                double fTimeSleeping = std::get<1>( vTimeThreadSleeping[ uiI ] ).count( );
                // if thread uiI is sleeping currently subtract the time it is currently sleeping from the awake time
                // since that amount will be added only once the thread wakes up...
                if( std::get<2>( vTimeThreadSleeping[ uiI ] ) )
                {
                    // explicit definition of xExtraSleepTime is required, so that .count() expresses the result in
                    // seconds.
                    duration xExtraSleepTime = xNow - std::get<0>( vTimeThreadSleeping[ uiI ] );
                    fTimeSleeping += xExtraSleepTime.count( );
                } // if
                xOut << uiI << "\t" << fTasksPerSecond << "\t" << vAmountBookedTasks[ uiI ] << "\t";

                xOut << 100 * ( 1.0 - fTimeSleeping / fTotalTime );
            } // if
            else
                // the queue that can be used by any thread requires a special print because for this queue
                // we do not have a sleep time entry as well as a amount booked entry
                xOut << "NO_T_ID\t" << fTasksPerSecond << "\tn/a\tn/a";
            // print the values for qlen=0 %qlen=1 %qlen=2 %qlen>2
            for( size_t uiJ = 0; uiJ < 4; uiJ++ )
                xOut << "\t" << 100 * vTimeAmountTasksInQueue[ uiI ][ uiJ ].count( ) / fTotalTime;
            xOut << std::endl;
        } // for
        // reset to default float output
        xOut << std::defaultfloat;
    } // method

    /**
     * @brief call this periodically to print a top like summary of the activity within the pool
     * @details
     * this calls printTime and then initTime every 10 seconds.
     */
    inline void top( std::ostream& xOut )
    {
        // return; // disable prints
        if( ( xLastTimePoint - xInitTimePoint ).count( ) > 30 /*amount of seconds between prints*/ )
        {
            printTime( xOut );
            initTime( vTimeThreadSleeping.size( ) );
        } // if
    } // method
#else
    // the following functions are empty dummies in case DO_PROFILE is set to zero

    inline void initTime( size_t uiPoolSize )
    {} // method
    inline void incTasks( size_t uiThreadId )
    {} // method
    inline void bookTask( size_t uiThreadId )
    {} // method
    inline void updateTime( )
    {} // method
    inline void goingToSleep( size_t uiThreadId )
    {} // method
    inline void wakingUp( size_t uiThreadId )
    {} // method
    inline void printTime( std::ostream& xOut )
    {} // method
    inline void top( std::ostream& xOut )
    {} // method
#endif

    std::mutex xQueueMutex; // mutex for control of concurrent access of the tasks queue
    std::condition_variable xCondition; // For wait, notify synchronization purposes.
    std::shared_ptr<std::mutex> pPoolLock; // pointer to mutex synchronizes concurrent activities in the pool
    bool bStop; // Flag for 'indicating end of pool operation'. If true, the pool operation is terminated.
    int iNextBookedThreadId = 0;

  public:
    const size_t uiPoolSize; // size of the pool

    SQLDBConPool( const SQLDBConPool& ) = delete; // no copies of thread pools

    /* TO DO: Configuration via JSON.
     * Notice: This constructor can throw exceptions. When an exception is thrown from a constructor, the object is not
     * considered instantiated, and therefore its destructor will not be called.
     */
    SQLDBConPool( size_t uiPoolSize = 1, // pass pool size here
                  const json& jDBConData = json{} ) // Connection Info as JSON; format depends on the backend
        : vqTasks( uiPoolSize + 1 ),
          pPoolLock( std::make_shared<std::mutex>( ) ), // create the pool lock
          bStop( false ), // stop must be false in the beginning
          uiPoolSize( uiPoolSize )
    {
        initTime( uiPoolSize );
        // Check size of pool
        if( uiPoolSize == 0 )
            throw std::runtime_error( "SQLDBConPool: The requested pool size must be greater than zero." );

        // Create all connection managers.
        // The creation of the connections has to be done sequentially with MySQL or connection creation can fail.
        for( size_t uiItr = 0; uiItr < uiPoolSize; ++uiItr )
            vConPool.emplace_back( std::make_shared<PooledSQLDBCon<DBImpl>>( pPoolLock, jDBConData ) );

        // Create an initialize all workers.
        for( size_t uiThreadId = 0; uiThreadId < uiPoolSize; ++uiThreadId )
            vWorkers.emplace_back( [this, uiThreadId] {
                // Treads are not allowed to throw exceptions or we face crashes.
                // Therefore the swallowing of exceptions.
                doNoExcept(
                    [&] {
                        // Get a reference to the connection manager belonging to this thread.
                        auto pConManager = vConPool[ uiThreadId ];
                        pConManager->uiThreadId = uiThreadId; // inform connector about corresponding taskid

                        // The loop continues until the stop flags indicates termination.
                        while( true )
                        {
                            // Synchronization of mutual access to the task queue.
                            // Start of scope of exclusive access for a single thread.
                            std::unique_lock<std::mutex> xLock( this->xQueueMutex );

                            // If there is no indication for termination and if there is nothing to do now,
                            // the lock is released so that some producer can push a fresh task into the queue.
                            // After inserting a fresh task, the producer calls notify_one() for awaking one of the
                            // waiting (sleeping) threads. This can be this thread or another in the pool.
                            while( !this->bStop && this->vqTasks[ this->uiPoolSize ].empty( ) &&
                                   this->vqTasks[ uiThreadId ].empty( ) )
                            {
                                goingToSleep( uiThreadId );
                                this->xCondition.wait( xLock );
                                // When we continue here, we have exclusive access once again.
                                wakingUp( uiThreadId );
                            } // while

                            // If there is an indication for termination and if there is nothing to do any more, we
                            // terminate the thread by returning.
                            if( this->bStop && this->vqTasks[ this->uiPoolSize ].empty( ) &&
                                this->vqTasks[ uiThreadId ].empty( ) )
                                return;

                            // We now have exclusive access to a non-empty task queue.
                            // Search for a matching task in the queues of tasks. Check the queue for this specific
                            // thread first, if there is no task in this queue then check the general queue.
                            size_t uiQueueId = uiThreadId;
                            if( this->vqTasks[ uiThreadId ].empty( ) )
                                uiQueueId = this->uiPoolSize;
                            assert( !this->vqTasks[ uiQueueId ].empty( ) );

                            updateTime( );

                            // We pick the matching task for execution and remove it from the respective queue.
                            TaskType fTaskToBeExecuted( this->vqTasks[ uiQueueId ].front( ) );
                            this->vqTasks[ uiQueueId ].pop( );

                            // Exclusive access to the task queue is not required any longer.
                            xLock.unlock( );

                            // Execute the task delivering its connection manager as argument.
                            // If the task throws an exception, this exception has to be caught via the future.
                            // (See the notice in the description of this class.)
                            fTaskToBeExecuted( pConManager );
                        } // while
                    }, // lambda ( doNoExcept )
                    // doNoExcept error message:
                    "SQLDBConPool:: The worker with Id " + std::to_string( uiThreadId ) +
                        "failed due to an internal error and terminated." );
            } ); // emplace_back
    } // constructor

    /** @brief Terminates the operation of the connection pool.
     *  Method blocks until all tasks in the queue have been executed.
     */
    void shutdown( )
    {
        {
            std::unique_lock<std::mutex> xLock( this->xQueueMutex );
            if( bStop )
                // Shutdown started already; no action required.
                // Return will release the mutex too ...
                return;
            bStop = true;
        } // end of scope for queue_mutex, so we unlock the mutex

        //std::cout << "Pool shutdown started ..." << std::endl;
        // Notify all waiting threads to continue so that they recognize the stop flag and terminate.
        xCondition.notify_all( );

        // Wait until all workers finished their job.
        // (The destructor is blocked until all workers finished their job.)
        // std::cout << "WORKER:" << vWorkers.size( ) << std::endl;
        for( size_t uiItr = 0; uiItr < vWorkers.size( ); ++uiItr )
            vWorkers[ uiItr ].join( );

        //std::cout << "Pool shutdown finished ..." << std::endl;
    } // method

    /** @brief The destructor blocks until all workers terminated. */
    ~SQLDBConPool( )
    {
        //std::cout << "SQLDBConPool Destructor." << std::endl;
        this->shutdown( );
        updateTime( );
        printTime( std::cout );
    } // destructor

    /** @brief Delivers a valid connection id for enqueuing for a specific worker. Using repeatedly the same worker is
     * reasonable, if two separate enqueues use the same query object.  */
    int getDedicatedConId( )
    {
        std::unique_lock<std::mutex> xLock( this->xQueueMutex );

        int iRet = iNextBookedThreadId;
        // increment and cycle (always leave at least one thread out of the cycle)
        iNextBookedThreadId = ( iNextBookedThreadId + 1 ) % ( uiPoolSize > 1 ? uiPoolSize - 1 : uiPoolSize );
        bookTask( iRet );
        return iRet;
    } // method

    const int NO_THREAD_ID = -1;

    /** @brief Enqueue a function func for execution via the DB pool using specific connection uiThreadId */
    template <class F, class... Args>
    auto inline enqueue( int iThreadId, F&& func, Args&&... args )
        -> std::future<typename std::result_of<F( std::shared_ptr<PooledSQLDBCon<DBImpl>>, Args... )>::type>
    {
        using return_type = typename std::result_of<F( std::shared_ptr<PooledSQLDBCon<DBImpl>>, Args... )>::type;

        // Don't allow enqueueing jobs after stopping the pool.
        {
            std::unique_lock<std::mutex> xLock( this->xQueueMutex );
            if( bStop )
                throw std::runtime_error( "SQLDBConPool: You tried to enqueue on stopped thread pool." );
        } // end of scope for queue_mutex, so we unlock the mutex

        // std::placeholders::_1 and represents some future value (this value will be only known
        // during function definition.
        // TO DO: Use an auto lambda of C++ 17
        auto xTask = std::make_shared<std::packaged_task<return_type( std::shared_ptr<PooledSQLDBCon<DBImpl>> )>> //
            ( std::bind( std::forward<F>( func ), std::placeholders::_1, std::forward<Args>( args )... ) );

        // The future outcome of the task. The caller will be later blocked until this result will
        // be available.
        std::future<return_type> xFuture = xTask->get_future( );
        {
            // Mutual access to the task queue has to be synchronized.
            std::unique_lock<std::mutex> lock( xQueueMutex );

            // Here we select the correct queue to put the task into.
            // There are two options: either a specific thread id is requested (iThreadId != NO_THREAD_ID)
            // then the task is put into the respective queue
            // or no specific thread is requested. then the task is put into the general queue
            if( iThreadId == NO_THREAD_ID )
                iThreadId = static_cast<int>( uiPoolSize );

            top( std::cout );
            updateTime( );
            incTasks( iThreadId );

            // The task gets as input a shared pointer to a DBConnection for doing its job
            vqTasks[ iThreadId ].push(
                [xTask]( std::shared_ptr<PooledSQLDBCon<DBImpl>> pDBCon ) { ( *xTask )( pDBCon ); } );
        } // end of scope of lock (lock gets released)

        if( iThreadId == static_cast<int>( uiPoolSize ) )
            // Inform one waiting workers that we have a fresh task. Because any worker can execute the task.
            this->xCondition.notify_one( );
        else
            // Inform all waiting workers that we have a fresh task. Only one worker will be able to execute the task
            this->xCondition.notify_all( );

        return xFuture;
    } // method enqueue

    /** @brief Enqueue a function func for execution via the DB pool */
    template <class F, class... Args>
    auto inline enqueue( F&& func, Args&&... args )
        -> std::future<typename std::result_of<F( std::shared_ptr<PooledSQLDBCon<DBImpl>>, Args... )>::type>
    {
        return enqueue( NO_THREAD_ID, std::forward<F>( func ), std::forward<Args>( args )... );
    } // method enqueue without thread id


    template <class F, class... Args>
    auto inline run( int iThreadId, F&& func, Args&&... args ) ->
        typename std::result_of<F( std::shared_ptr<PooledSQLDBCon<DBImpl>>, Args... )>::type
    {
        return enqueue( iThreadId, std::forward<F>( func ), std::forward<Args>( args )... ).get( );
    } // method

    template <class F, class... Args>
    auto inline run( F&& func, Args&&... args ) ->
        typename std::result_of<F( std::shared_ptr<PooledSQLDBCon<DBImpl>>, Args... )>::type
    {
        return run( NO_THREAD_ID, std::forward<F>( func ), std::forward<Args>( args )... );
    } // method

    // Usage example:
    //   dbpool.enqueue([&](std::shared_ptr<SQLDB<MySQLConDB>> pDBCon) { ... do something using the connection ... ;})
    // The statement delivers a future of type DBResponse that comprises the outcome of the operation.
    // Because threads are not allowed to terminate with still  a 'flying' exception, all exceptions must be catch in
    // the thread itself.
    // dbpool.
    //
}; // DB connection pool
