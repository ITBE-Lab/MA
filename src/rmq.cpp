#include "rmq.h"
#include <random>

void testRMQ()
{/*
    std::vector<RMQ<unsigned int>::RMQData> data = std::vector<RMQ<unsigned int>::RMQData>();
    for(unsigned int i = 1; i <= 1000; i++)
    {
        unsigned int x = std::rand() % 100;
        unsigned int y = std::rand() % 100;
        data.push_back(RMQ<unsigned int>::RMQData(x, y, i));
    }//for
    RMQ<unsigned int> rmq = RMQ<unsigned int>(data);

    for(int i = 1000; i > 0; i--)
    {
        assert(rmq.rmq(0,0, 100,100).score == i);
        rmq.rmq(0,0, 100,100).score = 0;
        std::cout << "check" << std::endl;
    }//for

    std::cout << ":D" << std::endl;
*/
}//function

void exportTest()
{
    boost::python::def("testRMQ", &testRMQ);
}//function