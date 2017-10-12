/** 
 * @file balancedSearchTree.h
 * @brief C++ Program to Implement Self Balancing Binary Search Tree
 * @note taken from: http://www.sanfoundry.com/cpp-program-implement-self-balancing-binary-tree/
 * @details
 *
 * The original code used pointers and did not care about free-ing anything.
 * Also only ints could be used as data before.
 * @note The data has to implement operator<= operator< and operator>.
 */

#ifndef BALANCED_SEARCH_TREE
#define BALANCED_SEARCH_TREE

#include <iostream>
#include <cstdlib>
#include <memory>

//forward declaration
template<class T>
class SelfBalancingBinarySearchTree;

/*
 * @brief Self balancing binary search tree node.
 */
template<class T>
class SBBSTNode
{
private:
    //allow the tree class to access all private stuff.
    friend class SelfBalancingBinarySearchTree<T>;
    int height;
    T* data;
    std::shared_ptr<SBBSTNode<T>> pLeft, pRight;
    std::weak_ptr<SBBSTNode<T>> pUp;

    SBBSTNode()
            :
        height(0),
        data(nullptr),
        pLeft(),
        pRight(),
        pUp()
    {}//constructor

    SBBSTNode(T* n)
            :
        height(0),
        data(n),
        pLeft(),
        pRight(),
        pUp()
    {}//constructor
public:
    std::shared_ptr<SBBSTNode<T>> next()
    {
        if(pRight == nullptr && pUp.expired())
            return nullptr;
        std::shared_ptr<SBBSTNode<T>> pRet = nullptr;
        if(pRight == nullptr)
        {
            SBBSTNode<T>* pDrag = this;
            pRet = pUp.lock();
            while(!pRet->pUp.expired() && (pRet->pRight == nullptr || pRet->pRight.get() == pDrag))
            {
                pDrag = pRet.get();
                pRet = pRet->pUp.lock();
            }//while
            if(pRet->pRight == nullptr || pRet->pRight.get() == pDrag)
                return nullptr;
            pRet = pRet->pRight;
        }//if
        else
            pRet = pRight;
        while(pRet != nullptr && pRet->pLeft != nullptr)
            pRet = pRet->pLeft;
        assert(pRet.get() != this);
        return pRet;
    }//function

    T* get()
    {
        return data;
    }//function
};//class

/*
 * @brief Self balancing binary search tree.
 */
template<class T>
class SelfBalancingBinarySearchTree
{
private:
    std::shared_ptr<SBBSTNode<T>> pRoot;



    /* Function to get height of node */
    int height(std::shared_ptr<SBBSTNode<T>> t) const
    {
        return t == nullptr ? -1 : t->height;
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

    /* Rotate binary tree node with pLeft child */
    std::shared_ptr<SBBSTNode<T>> rotateWithLeftChild(std::shared_ptr<SBBSTNode<T>> k2)
    {
        if(k2->pLeft == nullptr)
            return k2;
        std::shared_ptr<SBBSTNode<T>> k1 = k2->pLeft;
        k2->pLeft = k1->pRight;
        if(k1->pRight != nullptr)
            k1->pRight->pUp = k2;
        k1->pRight = k2;
        k1->pUp = k2->pUp;
        k2->pUp = k1;
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
        if(k2->pLeft != nullptr)
            k2->pLeft->pUp = k1;
        k2->pLeft = k1;
        k2->pUp = k1->pUp;
        k1->pUp = k2;
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
    std::shared_ptr<SBBSTNode<T>> insert(
            T* x,
            std::shared_ptr<SBBSTNode<T>> t,
            std::shared_ptr<SBBSTNode<T>>& pNewNode
        )
    {
        if (t == nullptr)
        {
            t = std::shared_ptr<SBBSTNode<T>>(new SBBSTNode<T>(x));
            pNewNode = t;
        }
        else if (*x <= *t->data)
        {
            DEBUG_3(
                std::cout << "\tinserting left" << std::endl;
            )
            t->pLeft = insert(x, t->pLeft, pNewNode);
            t->pLeft->pUp = t;
            if (height(t->pLeft) - height(t->pRight) == 2)
            {
                if (*x < *t->pLeft->data)
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
            t->pRight = insert(x, t->pRight, pNewNode);
            t->pRight->pUp = t;
            if (height(t->pRight) - height(t->pLeft) == 2)
            {
                if (*x > *t->pRight->data)
                    t = rotateWithRightChild(t);
                else
                    t = doubleWithRightChild(t);
            }
        }
        t->height = std::max(height(t->pLeft), height(t->pRight)) + 1;
        return t;
    }


public:
    SelfBalancingBinarySearchTree()
            :
        pRoot()
    {}//constructor

    /** @brief Function to check if tree is empty */
    bool isEmpty() const
    {
        return pRoot == nullptr;
    }

    /**
     * @brief Function to insert data.
     * @returns the next element of the tree.
     */
    std::shared_ptr<SBBSTNode<T>> insert(T* data)
    {
        std::shared_ptr<SBBSTNode<T>> pRet = nullptr;
        pRoot = insert(data, pRoot, pRet);
        assert(pRoot != nullptr);
        assert(pRet->pUp.lock() != pRet);
        return pRet;
    }

    /**
     * @brief Delete the leftmost element of the tree.
     */
    void deleteFirst()
    {
        pRoot = deleteFirst(pRoot);
    }

    /* @brief returns the leftmost element of the Tree. */
    T* first() const
    {
        std::shared_ptr<SBBSTNode<T>> r = pRoot;
        while (r->pLeft != nullptr)
                r = r->pLeft;
        return r->data;
    }
};//class

#endif