#include <iostream>
#include <vector>
#include <chrono>

#include "threadPool.h"
#include <boost/log/trivial.hpp>

void itemWorker( int i, size_t tid, int j )
{
	BOOST_LOG_TRIVIAL( trace ) << "Start data element" << i << " tid: " << tid << " j " << j;
	// std::cout << "hello " << i << std::endl;
	std::this_thread::sleep_for( std::chrono::seconds( 1 ) );
	BOOST_LOG_TRIVIAL( trace ) << "End data element"<< i << " tid: " << tid << " j " << j;
	// std::cout << "world " << i << std::endl;
	// return i*i;
}


int main_threadpool()
{
    std::vector<int> vDataItems = {5, 7, 8, 12, 45, 78, 5, 7, 90, 4, 2, 3, 1};

   
    std::vector< std::future<int> > results;
#if 0
    for(int i = 0; i < 8; ++i) {
        results.push_back(
            pool.enqueue([i] {
				BOOST_LOG_TRIVIAL(trace) << "Start data element" << i;
                // std::cout << "hello " << i << std::endl;
                std::this_thread::sleep_for(std::chrono::seconds(1));
				BOOST_LOG_TRIVIAL(trace) << "End data element" << i;
                // std::cout << "world " << i << std::endl;
                return i*i;
            })
        );
    }   
#endif
	// for (auto i : vDataItems )
	// {
	// 	results.push_back(
	// 		pool.enqueue( [i] ( size_t id, int j ){ return itemWorker(i, id, j); }, 6 ) );
	// }
	{
		ThreadPool pool( 4 );
		for ( int32_t i = 0; i < 10; i++ )
		{
			pool.enqueue( [i] ( size_t id, int j ) { itemWorker( i, id, j ); }, 6 );
		}
	} // Destruction of ThreadPool after all task terminated
	
	BOOST_LOG_TRIVIAL( trace ) << "*** CONTINUE ***";


	/* Because we work with futures we get the result here in this loop as soon as the corresponding task delivers it.
	 */
	for ( size_t i = 0; i < results.size(); ++i )
		BOOST_LOG_TRIVIAL( trace ) << "*** See result for " << i << " as " << results[i].get();
        // std::cout << results[i].get() << ' ';
    
	//for ( auto &value : results )
	//	BOOST_LOG_TRIVIAL( trace ) << "*** See result for " << value.get();

	std::cout << std::endl;
    
    return 0;
}
