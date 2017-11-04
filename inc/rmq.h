/**
 * @file rmq.h
 * @brief implements a n-dimensional range maximum query
 * @author Markus Schmidt
 */

#ifndef RMQ_H
#define RMQ_H

#include "searchTree.h"

class RMQData
{
public:
};//class


/**
 * @brief A Range Maximum Query implementation.
 * @details
 * A RMQ is a N dimensional binary search tree,
 * with canonical nodes as leafs in the N-th dimension.
 * @see Chaining algorithms for multiple genome comparison 
 * [Mohamed Ibrahim Abouelhoda, Enno Ohlebusch]
 */
class RMQ
{
private:
    class Node
    {
        std::vector<>
    };//class

    class Leaf : public Node
    {

    };//class

    class Branch : public Node
    {
        std::shared_ptr<Node> pL, pR;
    };//class


public:
};//class


#endif