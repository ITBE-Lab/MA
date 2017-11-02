/**
 * @file rmq.h
 * @brief implements a n-dimensional range maximum query
 * @author Markus Schmidt
 */

#ifndef RMQ_H
#define RMQ_H

#include "searchTree.h"

/**
 * @brief A Range Maximum Query implementation.
 * @details
 * A RMQ is a N dimensional binary search tree,
 * with canonical nodes as leafs in the N-th dimension.
 * @see Chaining algorithms for multiple genome comparison 
 * [Mohamed Ibrahim Abouelhoda, Enno Ohlebusch]
 */
template <unsigned int N>
class RMQ{
    SearchTree<RMQ<N-1>> tree;

};//class

template <>
class RMQ<0>{
    //SearchTree<> tree;

};//class

#endif