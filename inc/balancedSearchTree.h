#ifndef BALANCED_SEARCH_TREE
#define BALANCED_SEARCH_TREE
//teken from: http://www.sanfoundry.com/cpp-program-implement-self-balancing-binary-tree/
/*
 *  C++ Program to Implement Self Balancing Binary Search Tree
 */
#include <iostream>
#include <cstdlib>
#include <memory>

/* Class SBBSTNode */
template<class T>
class SBBSTNode
{
public:
    int height;
    T data;
    std::shared_ptr<SBBSTNode<T>> pLeft, pRight;

    SBBSTNode()
            :
        height(0),
        data(nullptr),
        pLeft(),
        pRight()
    {}//constructor

    SBBSTNode(T n)
            :
        height(0),
        data(n),
        pLeft(),
        pRight()
    {}//constructor
};//class

/*
 * Class SelfBalancingBinarySearchTree
 */

template<class T>
class SelfBalancingBinarySearchTree
{
private:
    std::shared_ptr<SBBSTNode<T>> pRoot;
public:
    SelfBalancingBinarySearchTree()
            :
        pRoot()
    {}//constructor

    /* Function to check if tree is empty */
    bool isEmpty() const
    {
        return pRoot == nullptr;
    }

        /* Function to get height of node */
    int height(std::shared_ptr<SBBSTNode<T>> t) const
    {
        return t == nullptr ? -1 : t->height;
    }

    /* Rotate binary tree node with pLeft child */
    std::shared_ptr<SBBSTNode<T>> rotateWithLeftChild(std::shared_ptr<SBBSTNode<T>> k2)
    {
        if(k2->pLeft == nullptr)
            return k2;
        std::shared_ptr<SBBSTNode<T>> k1 = k2->pLeft;
        k2->pLeft = k1->pRight;
        k1->pRight = k2;
        k2->height = std::max(height(k2->pLeft), height(k2->pRight)) + 1;
        k1->height = std::max(height(k1->pLeft), k2->height) + 1;
        return k1;
    }

    /* Rotate binary tree node with pRight child */
    std::shared_ptr<SBBSTNode<T>> rotateWithRightChild(std::shared_ptr<SBBSTNode<T>> k1)
    {
        if(k1->pRight == nullptr)
            return k1;
        std::shared_ptr<SBBSTNode<T>> k2 = k1->pRight;
        k1->pRight = k2->pLeft;
        k2->pLeft = k1;
        k1->height = std::max(height(k1->pLeft), height(k1->pRight)) + 1;
        k2->height = std::max(height(k2->pRight), k1->height) + 1;
        return k2;
    }

    /*
    * Double rotate binary tree node: first pLeft child
    * with its pRight child; then node k3 with new pLeft child
    */
    std::shared_ptr<SBBSTNode<T>> doubleWithLeftChild(std::shared_ptr<SBBSTNode<T>> k3)
    {
        k3->pLeft = rotateWithRightChild(k3->pLeft);
        return rotateWithLeftChild(k3);
    }

    /*
    * Double rotate binary tree node: first pRight child
    * with its pLeft child; then node k1 with new pRight child
    */
    std::shared_ptr<SBBSTNode<T>> doubleWithRightChild(std::shared_ptr<SBBSTNode<T>> k1)
    {
        k1->pRight = rotateWithLeftChild(k1->pRight);
        return rotateWithRightChild(k1);
    }

    /* Function to insert data recursively */
    std::shared_ptr<SBBSTNode<T>> insert(T x, std::shared_ptr<SBBSTNode<T>> t, T*& pNext)
    {
        if (t == nullptr)
            t = std::shared_ptr<SBBSTNode<T>>(new SBBSTNode<T>(x));
        else if (x <= t->data)
        {
            DEBUG_3(
                std::cout << "\tinserting left" << std::endl;
            )
            t->pLeft = insert(x, t->pLeft, pNext);
            if(pNext == nullptr)
                pNext = &t->data;
            if (height(t->pLeft) - height(t->pRight) == 2)
            {
                if (x < t->pLeft->data)
                    t = rotateWithLeftChild(t);
                else
                    t = doubleWithLeftChild(t);
            }
        }
        else
        {
            DEBUG_3(
                std::cout << "\tinserting right" << std::endl;
            )
            t->pRight = insert(x, t->pRight, pNext);
            if (height(t->pRight) - height(t->pLeft) == 2)
            {
                if (x > t->pRight->data)
                    t = rotateWithRightChild(t);
                else
                    t = doubleWithRightChild(t);
            }
        }
        t->height = std::max(height(t->pLeft), height(t->pRight)) + 1;
        return t;
    }

    /* Function to insert data 
     * returns the next element of the tree
     */
    T* insert(T data)
    {
        T* pRet = nullptr;
        pRoot = insert(data, pRoot, pRet);
        assert(pRoot != nullptr);
        return pRet;
    }

    
    /* Function to insert data recursively */
    std::shared_ptr<SBBSTNode<T>> deleteFirst(std::shared_ptr<SBBSTNode<T>> t)
    {
        if (t == nullptr)
            return nullptr;
        else if (t->pLeft == nullptr)
            return t->pRight;
        else
        {
            //TODO: no rotation needed?
            t->pLeft = deleteFirst(t->pLeft);
        }
        t->height = std::max(height(t->pLeft), height(t->pRight)) + 1;
        return t;
    }

    /* Function to insert data 
     * returns the next element of the tree
     */
    void deleteFirst()
    {
        pRoot = deleteFirst(pRoot);
    }

    /* Functions to search for an element */
    T& first() const
    {
        std::shared_ptr<SBBSTNode<T>> r = pRoot;
        while (r->pLeft != nullptr)
                r = r->pLeft;
        return r->data;
    }
#if 0
    T search(std::shared_ptr<SBBSTNode<T>> r, int val)
    {
        while (r != nullptr)
        {
            if (r->data > val)
                r = r->pLeft;
            else if (r->data < val)
                r = r->pRight;
            else
                return r->data;
        }//while
        return nullptr;
    }

    /* Functions to search for an element */
    T search(int val)
    {
        return search(pRoot, val);
    }
#endif
};//class

#endif