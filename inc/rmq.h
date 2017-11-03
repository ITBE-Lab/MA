/**
 * @file rmq.h
 * @brief implements a n-dimensional range maximum query
 * @author Markus Schmidt
 */

#ifndef RMQ_H
#define RMQ_H

#include "searchTree.h"

class RMQData{
public:
    virtual bool smaller(const RMQData& other, unsigned int iDimension) const;
    virtual bool equal(const RMQData& other, unsigned int iDimension) const;

    bool operator<(const RMQData& other) const
    {
        return smaller(other, 0);
    }//function
    bool operator==(const RMQData& other) const
    {
        return equal(other, 0);
    }//function

    virtual unsigned int numDimensions() const;
};//class

bool operator<(const std::shared_ptr<RMQData> a, const std::shared_ptr<RMQData> b) const
{
    return *a < *b;
}//function

bool operator==(const std::shared_ptr<RMQData> a, const std::shared_ptr<RMQData> b) const
{
    return *a == *b;
}//function


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
    SearchTree<std::shared_ptr<RMQ<N-1>>> tree;

public:
    RMQ(std::shared_ptr<RMQData> data)
            :
        tree()
    {
        tree.insert(std::shared_ptr<RMQ<N-1>>(new RMA<N-1>(data)));
    }//constructor

    RMQ()
            :
        tree()
    {}//constructor

    const RMQData& root() const
    {
        //get the root of the tree which is itself a RMQ; call root() on that.
        return tree.root().root();
    }//function

    bool operator<(const RMQ<N>& other) const
    {
        return root().smaller(other.root(), N);
    }//function

    bool operator==(const RMQ<N>& other) const
    {
        return root().equal(other.root(), N);
    }//function

    void insert(std::shared_ptr<RMQData> data)
    {
        //sanity check
        assert(data->numDimensions() >= N);

        std::shared_ptr<RMQ<N-1>> newSubTree = std::shared_ptr<RMQ<N-1>>(new RMQ<N-1>(data));
        //returns an iterator thus the * in order to dereference it
        std::shared_ptr<RMQ<N-1>> foundTree = *tree.insertOrGet(newSubTree);

        //insert the data into the found subtree
        if(foundTree != newSubTree)
        {
            foundTree->insert(data);
        }//if
        //the else case is done by creating the newSubTree with the data already inside.

    }//function

};//class

template <unsigned int N>
bool operator<(const std::shared_ptr<RMQ<N>> a, const std::shared_ptr<RMQ<N>> b) const
{
    return *a < *b;
}//function

template <unsigned int N>
bool operator ==(const std::shared_ptr<RMQ<N>> a, const std::shared_ptr<RMQ<N>> b) const
{
    return *a == *b;
}//function

template <>
class RMQ<1>{
    SearchTree<std::shared_ptr<RMQData>> tree;

public:
    RMQ(std::shared_ptr<RMQData> data)
            :
        tree()
    {
        tree.insert(data);
    }//constructor

    const RMQData& root() const
    {
        return tree.root();
    }//function

    void insert(std::shared_ptr<RMQData> data)
    {
        tree.insert(data);
    }//function

};//class

#endif