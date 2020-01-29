/* Authors: Arne Kutzner and Markus Schmidt
 * Created: Jan. 2020
 * This file is part of the ITBE-Lab code collection.
 * MIT License
 * @file db_con_pool.h
 * @brief Connection pool for databases
 */

#pragma once

#include <condition_variable>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

#include "common.h" // general API for SQL like DB

/* Basic characteristics:
 * - The pool synchronizes threads that ask for access, but the pool itself does not comprises threads.
 *   (Therefore, e.g. an additional thread pool is required for delivering multiple thread)
 */
template <typename DBImpl> class SQLDBConPool
{
  public:
    class ConnectionManager
    {
      private:
        size_t uiTaskId; // id of the task belonging to the ConnectionManager

      public:
        ConnectionManager( const ConnectionManager& ) = delete; // no copies of connection managers
        std::shared_ptr<SQLDB<DBImpl>> pDBCon;

        /** @brief Constructs a DB connection manager.
         */
        ConnectionManager( ) : pDBCon( std::make_shared<SQLDB<DBImpl>>( ) )
        {
            std::cout << "Create ConnectionManager " << std::endl;
        } // constructor

        inline size_t getTaskId( )
        {
            return this->uiTaskId;
        } // method

        friend class SQLDBConPool<DBImpl>;
    }; // class

  private:
    using TaskType =
        std::function<void( std::shared_ptr<ConnectionManager> )>; // type of the tasks to be executed

    std::vector<std::thread> vWorkers; // worker threads; each thread manages one connection
    std::vector<std::shared_ptr<ConnectionManager>> vConPool; // pointers to actual connections
    std::queue<TaskType> qTasks; // queue of lambda functions that shall to executed
    std::mutex xQueueMutex; // mutex for control of concurrent access of the tasks queue
    std::condition_variable xCondition; // For wait, notify synchronization purposes.
    bool bStop; // Flag for 'indicating end of pool operation'. If true, the pool operation is terminated.

  public:
    const size_t uiPoolSize; // size of the pool

    SQLDBConPool( const SQLDBConPool& ) = delete; // no copies of thread pools

    /* TO DO: Configuration via JSON.
     */
    SQLDBConPool( size_t uiPoolSize = 1 )
        : bStop( false ), // stop must be false in the beginning
          uiPoolSize( uiPoolSize )
    {
        // Check size of pool
        if( uiPoolSize == 0 )
            throw std::runtime_error( "SQLDBConPool: The requested pool size must be greater than zero." );

        // Create all connection managers.
        // The creation of the connections has to be done sequentially with MySQL or connection creation can fail.
        for( size_t uiItr = 0; uiItr < uiPoolSize; ++uiItr )
            vConPool.emplace_back( std::make_shared<ConnectionManager>( ) );

        // Create an initialize all workers.
        for( size_t uiTaskId = 0; uiTaskId < uiPoolSize; ++uiTaskId )
            vWorkers.emplace_back( [ this, uiTaskId ] {
                // Treads in a pool are not allowed to throw exceptions or we face crashes.
                // Therefore the swallowing of exceptions.
                doSwallowingExcpt( [ & ] {
                    // Get a reference to the connection manager belonging to this thread.
                    auto pConManager = vConPool[ uiTaskId ];
                    pConManager->uiTaskId = uiTaskId; // inform connector about corresponding taskid

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
                        while( !this->bStop && this->qTasks.empty( ) )
                        {
                            this->xCondition.wait( xLock );
                            // When we continue here, we have exclusive access once again.
                        } // while

                        // If there is an indication for termination and if there is nothing to do any more, we
                        // terminate the thread by returning.
                        if( this->bStop && this->qTasks.empty( ) )
                        {
                            std::cout << "Thread " << uiTaskId << " terminated" << std::endl;
                            return;
                        } // if

                        // We now have exclusive access to a non-empty task queue.
                        // We pick the front element (front lambda expression) for execution and remove it from the
                        // queue.
                        TaskType fTaskToBeExcuted( this->qTasks.front( ) );
                        this->qTasks.pop( );

                        // Exclusive access to the task queue is not required any longer.
                        xLock.unlock( );

                        // Execute the task delivering its connection manager as argument.
                        // Swallow exceptions so that the thread can continue even if something goes wrong.
						// MARKUS: Encapsualte into try-catch ....
                        doSwallowingExcpt( [ & ] (){ fTaskToBeExcuted( pConManager ); } );
                    } // while

                    // If we ever come to this point, a serious problem occurred.
                } ); // lambda
            } ); // emplace_back
    } // constructor

    /** brief@ The destructor blocks until all workers terminated. */
    ~SQLDBConPool( )
    {
        std::cout << "Pool termination start" << std::endl;
        // Indicate end of pool operation by setting bStop to true.
        {
            std::unique_lock<std::mutex> xLock( this->xQueueMutex );
            bStop = true;
        } // end of scope for queue_mutex, so we unlock the mutex

        // Notify all waiting threads to continue so that they recognize the stop flag and terminate.
        xCondition.notify_all( );

        // Wait until all workers finished their job.
        // (The destructor is blocked until all workers finished their job.)
        for( size_t uiItr = 0; uiItr < vWorkers.size( ); ++uiItr )
            vWorkers[ uiItr ].join( );
        std::cout << "Pool termination end" << std::endl;
    } // destructor

    /** @brief Enqueue a function func for execution via the DB pool */
    template <class F, class... Args>
    auto inline enqueue( F&& func, Args&&... args )
        -> std::future<typename std::result_of<F(std::shared_ptr<ConnectionManager>, Args... )>::type>
    {
        using return_type = typename std::result_of<F(std::shared_ptr<ConnectionManager>, Args... )>::type;

        // Don't allow enqueueing jobs after stopping the pool.
        if( bStop )
            throw std::runtime_error( "SQLDBConPool: You tried to enqueue on stopped thread pool." );

        // std::placeholders::_1 and represents some future value (this value will be only known
        // during function definition.
        // TO DO: Use an auto lambda of C++ 17
        auto xTask = std::make_shared<std::packaged_task<return_type(std::shared_ptr<ConnectionManager>)>> //
            ( std::bind( std::forward<F>( func ), std::placeholders::_1, std::forward<Args>( args )... ) );

        // The future outcome of the task. The caller will be later blocked until this result will
        // be available.
        std::future<return_type> xFuture = xTask->get_future( );
        {
            // Mutual access to the task queue has to be synchronized.
            std::unique_lock<std::mutex> lock( xQueueMutex );

            // The task gets as input a shared pointer to a DBConnection for doing its job
            qTasks.push( [ xTask ](std::shared_ptr<ConnectionManager> pDBCon ) { ( *xTask )( pDBCon ); } );
        } // end of scope of lock (lock gets released)
        // Inform some waiting consumer (worker) that we have a fresh task.
        this->xCondition.notify_one( );

        return xFuture;
    } // method enqueue

    /** @brief Executes the func so that exceptions are swallowed.
     *  @detail For destructors and threads so that they do not throw.
     */
    template <typename F> inline void doSwallowingExcpt( F&& func )
    {
        try
        {
			func( );
        } // try
        catch( std::exception& rxExcpt )
        {
            // Swallow the exception
            std::cout << std::string( "SQLDBConPool: Task/thread execution terminated abnormally. Details:\n" ) +
                             rxExcpt.what( )
                      << std::endl;
        } // catch
        catch( ... )
        {
            std::cout << "SQLDBConPool: Task/thread execution terminated abnormally for unknown reasons." << std::endl;
        } // catch
		std::cout << "Do swallowing terminated..." << std::endl;
    } // method

    // Usage example:
    //   dbpool.enqueue([&](std::shared_ptr<SQLDB<MySQLConDB>> pDBCon) { ... do something using the connection ... ;})
    // The statement delivers a future of type DBResponse that comprises the outcome of the operation.
    // Because threads are not allowed to terminate with still  a 'flying' exception, all exceptions must be catch in
    // the thread itself.
    // dbpool.
    //
}; // DB connection pool
