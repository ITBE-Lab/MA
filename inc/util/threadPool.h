/** 
 * @file threadPool.h
 * @brief A thread pool.
 * @author Jakob Progsch, Václav Zeman, Arne Kutzner, Markus Schmidt
 * @copyright
 * Copyright (c) 2012 Jakob Progsch, Václav Zeman
 *
 *
 * This software is provided 'as-is', without any express or implied
 * warranty. In no event will the authors be held liable for any damages
 * arising from the use of this software.
 * 
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 * 
 * 1. The origin of this software must not be misrepresented; you must not
 * claim that you wrote the original software. If you use this software
 * in a product, an acknowledgment in the product documentation would be
 * appreciated but is not required.
 * 
 * 2. Altered source versions must be plainly marked as such, and must not be
 * misrepresented as being the original software.
 * 
 * 3. This notice may not be removed or altered from any source
 * distribution.
 *
 * @details
 * Code originates from:
 * https://github.com/progschj/ThreadPool
 * Changes:
 * - The thread id extension was added by Arne Kutzner.
 * - ThreadPoolAllowingRecursiveEnqueue was added by Markus Schmidt.
 */


/* ORIGINAL SOURCE FROM https://github.com/progschj/ThreadPool :

    #ifndef THREAD_POOL_H
    #define THREAD_POOL_H

    #include <vector>
    #include <queue>
    #include <memory>
    #include <thread>
    #include <mutex>
    #include <condition_variable>
    #include <future>
    #include <functional>
    #include <stdexcept>

    class ThreadPool {
    public:
        ThreadPool(size_t);
        template<class F, class... Args>
        auto enqueue(F&& f, Args&&... args) 
            -> std::future<typename std::result_of<F(Args...)>::type>;
        ~ThreadPool();
    private:
        // need to keep track of threads so we can join them
        std::vector< std::thread > workers;
        // the task queue
        std::queue< std::function<void()> > tasks;
        
        // synchronization
        std::mutex queue_mutex;
        std::condition_variable condition;
        bool stop;
    };
    
    // the constructor just launches some amount of workers
    inline ThreadPool::ThreadPool(size_t threads)
        :   stop(false)
    {
        for(size_t i = 0;i<threads;++i)
            workers.emplace_back(
                [this]
                {
                    for(;;)
                    {
                        std::function<void()> task;

                        {
                            std::unique_lock<std::mutex> lock(this->queue_mutex);
                            this->condition.wait(lock,
                                [this]{ return this->stop || !this->tasks.empty(); });
                            if(this->stop && this->tasks.empty())
                                return;
                            task = std::move(this->tasks.front());
                            this->tasks.pop();
                        }

                        task();
                    }
                }
            );
    }

    // add new work item to the pool
    template<class F, class... Args>
    auto ThreadPool::enqueue(F&& f, Args&&... args) 
        -> std::future<typename std::result_of<F(Args...)>::type>
    {
        using return_type = typename std::result_of<F(Args...)>::type;

        auto task = std::make_shared< std::packaged_task<return_type()> >(
                std::bind(std::forward<F>(f), std::forward<Args>(args)...)
            );
            
        std::future<return_type> res = task->get_future();
        {
            std::unique_lock<std::mutex> lock(queue_mutex);

            // don't allow enqueueing after stopping the pool
            if(stop)
                throw std::runtime_error("enqueue on stopped ThreadPool");

            tasks.emplace([task](){ (*task)(); });
        }
        condition.notify_one();
        return res;
    }

    // the destructor joins all threads
    inline ThreadPool::~ThreadPool()
    {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            stop = true;
        }
        condition.notify_all();
        for(std::thread &worker: workers)
            worker.join();
    }

    #endif
*/

#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#define NUM_THREADS_ALIGNER 0



/*
 * IDEA: 
 * if all within module synchronization is in "SYNC()" blocs
 * i can enable and disable inner module multithreading simply by setting NUM_THREADS_ALIGNER to 0
 * 
 * in that case threadpool is set up to execute all tasks using the main thread.
 * WARNING:
 * the multithreading within the pledge class is still active
 * so no SYNC blocks in the threads class
 */
#if NUM_THREADS_ALIGNER > 0
#define SYNC(x) x
#else //DEBUG_LEVEL
#define SYNC(x)
#endif //DEBUG_LEVEL

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>

/**
 * @brief A Threadpool.
 * @details
 * // create thread pool with 4 worker threads\n
 * ThreadPool pool(4);\n
 *
 * // enqueue and store future\n
 * auto result = pool.enqueue( [](int answer) { return answer; }, 42 );\n
 *
 * // get result from future\n
 * std::cout << result.get() << std::endl;\n
 *
 *
 * if the threadpool is set to have 0 threads, we will execute every task in the main thread
 * immediately, thus we dont need to setup threads at all.
 * this can be used to disable multithreading simply setting the pool to have 0 threads
 */
class ThreadPool {
 private:
    /* need to keep track of threads so we can join them
	 */
    std::vector< std::thread > workers;
    
	/* the task queue
	 */
    std::queue< std::function<void(size_t)> > tasks;
    
    /* Synchronization
	 */
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool bStop;
	const size_t threads;

public:

	/* External definition
	 */
    ThreadPool( size_t );

	/* External definition
	 */
	~ThreadPool();

	/* add new work item to the pool
	 */
	template<class F, class... Args>
	auto enqueue( F &&f, Args&&... args ) 
		->  std::future<typename std::result_of<F( size_t, Args... )>::type>
	{
		typedef typename std::result_of<F( size_t, Args... )>::type return_type;

		/* Don't allow enqueueing after stopping the pool 
		 */
		if ( bStop )
		{
			throw std::runtime_error( "enqueue on stopped ThreadPool" );
		} // if

		/* std::placeholders::_1 and represents some future value (this value will be only known during function definition.
		 */
		auto task = std::make_shared< std::packaged_task< return_type( size_t ) > >
					(
						std::bind( std::forward<F>( f ), std::placeholders::_1, std::forward<Args>( args )... ) 
					);

		/* The future outcome of the task. The caller will be later blocked until this result will be available. 
		 */
		std::future<return_type> xFuture = task->get_future();

		//if the threadpool is set to have 0 threads then execute the task in the main thread
		if(threads == 0)
		{
			(*task)(0);
			return xFuture;
		}//if
		
		{
			/* Mutual access to the task queue has to be synchronized. 
			 */
			std::unique_lock<std::mutex> lock( queue_mutex );

			/* The task will get delivered its task_id later during its execution. 
			 */
			tasks.push( [task] (size_t task_id) { (*task)( task_id ); } );
		} // end of scope for the lock, so we will release the lock.

		/* Inform some waiting consumer (worker )that we have a fresh task.
		 */
		condition.notify_one();
		
		return xFuture;
	} // method enqueue
   
};
 
/* Constructor just launches some amount of workers
 */
inline ThreadPool::ThreadPool( size_t threads )
    	:
	bStop( false ), // stop must be false in the beginning
	threads(threads)
{
	//if the threadpool is set to have 0 threads, we will execute every task in the main thread
	//immediately, thus we dont need to setup threads at all.
	//this can be used to disable multithreading simply setting the pool to have 0 threads
	if(threads == 0)
		return;
	for ( size_t i = 0; i < threads; ++i )
        workers.emplace_back(
            [this, i]
            {
                /* So long we have some data the thread processes these data 
				 */
				for(;;)
                {
					/* Synchronization of mutual access to the task queue.
					 */
					std::unique_lock<std::mutex> lock( this->queue_mutex );

					while ( !this->bStop && this->tasks.empty() )
					{
						/* We release the lock, so that some producer can push some task into the queue.
						 * The produser will call notify_one(), in order to release 
						 */
						this->condition.wait( lock );
					} // while

					if ( this->bStop && this->tasks.empty() )
					{
						/* All work done and signal for termination received.
						 */
						return;
					} // if
					
					/* Initialize the function variable task, so that it refers the task at the top of the queue.
					 */
					std::function<void( size_t )> task( this->tasks.front() );
					this->tasks.pop();
					
                    lock.unlock();

					try
					{
						/* Execute the task (that we received as top of the queue)
						*/
						task( i );
					} 
					catch(std::exception e) 
					{
						std::cerr << "exception when executing task:" << std::endl;
						std::cerr << e.what() << std::endl;
					}
					catch(...) 
					{
						std::cerr << "unknown exception when executing task" << std::endl;
					}
                } // for
            } // lambda
        ); // function call
} // method

/* the destructor joins all threads
 */
inline ThreadPool::~ThreadPool()
{
	//if the threadpool is set to have 0 threads, we will execute every task in the main thread
	//immediately, thus we dont need to setup threads at all.
	//this can be used to disable multithreading simply setting the pool to have 0 threads
	if(threads == 0)
		return;
    {
        std::unique_lock<std::mutex> lock( queue_mutex );
        bStop = true;
    } // end of scope for queue_mutex, so we unlock the mutex
    condition.notify_all();

	/* We wait until all workers finished their job.
	 * (Destruction thread is blocked until all workers finished their job.)
	 */
	for ( size_t i = 0; i < workers.size(); ++i )
	{
		workers[i].join();
	} // for
}


/**
 * @brief A Threadpool.
 * @details
 * // create thread pool with 4 worker threads\n
 * ThreadPool pool(4);\n
 *
 * // enqueue and store future\n
 * auto result = pool.enqueue( [](int answer) { return answer; }, 42 );\n
 *
 * // get result from future\n
 * std::cout << result.get() << std::endl;\n
 *
 *
 * if the threadpool is set to have 0 threads, we will execute every task in the main thread
 * immediately, thus we dont need to setup threads at all.
 * this can be used to disable multithreading simply setting the pool to have 0 threads
 *
 * @note This pool allows enqueues from within a worker thread.
 */
class ThreadPoolAllowingRecursiveEnqueue {

private:
	/* need to keep track of threads so we can join them
	*/
	std::vector< std::thread > workers;

	/* the task queue
	*/
	std::queue< std::function<void(size_t)> > tasks;

	/* Synchronization
	*/
	std::mutex queue_mutex;
	std::condition_variable condition;
	bool bStop;
	const size_t threads;

public:

	/* External definition
	*/
	ThreadPoolAllowingRecursiveEnqueue(size_t);

	/* External definition
	*/
	~ThreadPoolAllowingRecursiveEnqueue();

	/* add new work item to the pool
	*/
	template<class F, class... Args>
	auto enqueue(F &&f, Args&&... args)
		->  std::future<typename std::result_of<F(size_t, Args...)>::type>
	{
		typedef typename std::result_of<F(size_t, Args...)>::type return_type;

		/* 
		 * std::placeholders::_1 and represents some future value 
		 * (this value will be only known during function definition.
		 */
		auto task = std::make_shared< std::packaged_task< return_type(size_t) > >
			(
			std::bind(std::forward<F>(f), std::placeholders::_1, std::forward<Args>(args)...)
			);

		/* The future outcome of the task. The caller will be later blocked until this result will be available.
		*/
		std::future<return_type> xFuture = task->get_future();

		//if the threadpool is set to have 0 threads then execute the task in the main thread
		if(threads == 0)
		{
			(*task)(0);
			return xFuture;
		}//if

		{
			/* Mutual access to the task queue has to be synchronized.
			*/
			std::unique_lock<std::mutex> lock(queue_mutex);

			/* The task will get delivered its task_id later during its execution.
			*/
			tasks.push([task](size_t task_id) { (*task)(task_id); });
		} // end of scope for the lock, so we will release the lock.

		/* Inform some waiting consumer (worker )that we have a fresh task.
		*/
		condition.notify_one();

		return xFuture;
	} // method enqueue
};

/* Constructor just launches some amount of workers
*/
inline ThreadPoolAllowingRecursiveEnqueue::ThreadPoolAllowingRecursiveEnqueue(size_t threads)
		:
	bStop( false ), // stop must be false in the beginning
	threads(threads)
{
	//if the threadpool is set to have 0 threads, we will execute every task in the main thread
	//immediately, thus we dont need to setup threads at all.
	//this can be used to disable multithreading simply setting the pool to have 0 threads
	if(threads == 0)
		return;
	for (size_t i = 0; i < threads; ++i)
		workers.emplace_back(
		[this, i]
	{
		/* So long we have some data the thread processes these data
		*/
		for (;;)
		{
			/* Synchronization of mutual access to the task queue.
			*/
			std::unique_lock<std::mutex> lock(this->queue_mutex);

			while (!this->bStop && this->tasks.empty())
			{
				/* We release the lock, so that some producer can push some task into the queue.
				* The produser will call notify_one(), in order to release
				*/
				this->condition.wait(lock);
			} // while

			if (this->bStop && this->tasks.empty())
			{
				/* All work done and signal for termination received.
				*/
				return;
			} // if

			/* Initialize the function variable task, so that it refers the task at the top of the queue.
			*/
			std::function<void(size_t)> task(this->tasks.front());
			this->tasks.pop();

			lock.unlock();

			/* Execute the task (that we received as top of the queue)
			*/
			task(i);
		} // for
	} // lambda
	); // function call
} // method

/* the destructor joins all threads
*/
inline ThreadPoolAllowingRecursiveEnqueue::~ThreadPoolAllowingRecursiveEnqueue()
{
	//if the threadpool is set to have 0 threads, we will execute every task in the main thread
	//immediately, thus we dont need to setup threads at all.
	//this can be used to disable multithreading simply setting the pool to have 0 threads
	if(threads == 0)
		return;

	{
		std::unique_lock<std::mutex> lock(queue_mutex);
		bStop = true;
	} // end of scope for queue_mutex, so we unlock the mutex
	condition.notify_all();

	/* We wait until all workers finished their job.
	* (Destruction thread is blocked until all workers finished their job.)
	*/
	for (size_t i = 0; i < workers.size(); ++i)
	{
		workers[i].join();
	} // for
}


#endif
