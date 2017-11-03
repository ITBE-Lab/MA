/** 
 * @file searchTree.h
 * @brief A Self Balancing Binary Search Tree.
 */

#ifndef BALANCED_SEARCH_TREE
#define BALANCED_SEARCH_TREE

#include <iostream>
#include <cstdlib>
#include <memory>
#include "exception.h"


/**
 * @brief Self balancing binary search tree.
 */
template<typename T>
class SearchTree
{
private:
    class Leaf;
    class Branch;

    /**
     * @brief Self balancing binary search tree node.
     */
    class Node
    {
    public:
        virtual unsigned int getHeight() const
        {
            throw NullPointerException("there should never be an instance of the node class");
        }
    private:
        std::shared_ptr<Node> pLeft, pRight;
        std::weak_ptr<Branch> pPrev, pNext;

        virtual void setHeight(unsigned int h)
        {
            throw NullPointerException("there should never be an instance of the node class");
        }
        void reCalcHeight()
        {
            setHeight(std::max(pLeft->getHeight(), pRight->getHeight())+1);
        }//function

        virtual bool operator<(T& data) const
        {
            throw NullPointerException("there should never be an instance of the node class");
        }

    public:
        Node() = delete;
        Node(
                std::shared_ptr<SearchTree<T>::Branch> pPrev,
                std::shared_ptr<SearchTree<T>::Branch> pNext,
                std::shared_ptr<SearchTree<T>::Node> pLeft,
                std::shared_ptr<SearchTree<T>::Node> pRight
            )
                :
            pLeft(pLeft),
            pRight(pRight),
            pPrev(pPrev),
            pNext(pNext)
        {}//constructor

        virtual const T& get() const
        {
            throw NullPointerException("there should never be an instance of the node class");
        }


        virtual std::shared_ptr<Branch> next() const
        {
            return pNext.lock();
        }//function

        virtual std::shared_ptr<Node> insert(
                T& data, 
                std::shared_ptr<Node> pThis,
                std::shared_ptr<Branch> pLastPrev,
                std::shared_ptr<Branch> pLastNext,
                std::shared_ptr<Branch>& pNew
            )
        {
            std::shared_ptr<Branch> pThis_c = std::dynamic_pointer_cast<Branch>(pThis);
            if(*pThis < data)
                pRight = pRight->insert(data, pRight, pThis_c, pLastNext, pNew);
            else
                pLeft = pLeft->insert(data, pLeft, pLastPrev, pThis_c, pNew);
            reCalcHeight();
            return pThis;
        }//function

        virtual std::shared_ptr<Node> insertOrGet(
                T& data, 
                std::shared_ptr<Node> pThis,
                std::shared_ptr<Branch> pLastPrev,
                std::shared_ptr<Branch> pLastNext,
                std::shared_ptr<Branch>& pNew
            )
        {
            std::shared_ptr<Branch> pThis_c = std::dynamic_pointer_cast<Branch>(pThis);
            if(*pThis == data)
            {
                pNew = pThis_c;
                return pThis;
            }//if
            else if(*pThis < data)
                pRight = pRight->insert(data, pRight, pThis_c, pLastNext, pNew);
            else
                pLeft = pLeft->insert(data, pLeft, pLastPrev, pThis_c, pNew);
            reCalcHeight();
            return pThis;
        }//function

        virtual std::shared_ptr<Node> rotateWithRightChild(std::shared_ptr<Node> pThis)
        {
            std::shared_ptr<Node> pSwapWith = getRightOrLeaf(pThis);
            pRight = pSwapWith->getLeftOrLeaf(pSwapWith);
            pSwapWith->pLeft = pThis;

            pRight->reCalcHeight();
            reCalcHeight();
            pSwapWith->reCalcHeight();

            return pSwapWith;
        }//function
        
        virtual std::shared_ptr<Node> rotateWithLeftChild(std::shared_ptr<Node> pThis)
        {
            std::shared_ptr<Node> pSwapWith = getLeftOrLeaf(pThis);
            pLeft = pSwapWith->getRightOrLeaf(pSwapWith);
            pSwapWith->pRight = pThis;

            pLeft->reCalcHeight();
            reCalcHeight();
            pSwapWith->reCalcHeight();
            return pSwapWith;
        }//function

        virtual std::shared_ptr<Node> doubleWithRightChild(std::shared_ptr<Node> pThis)
        {
            pLeft = pLeft->rotateWithRightChild(pLeft);
            return rotateWithLeftChild(pThis);
        }//function

        virtual std::shared_ptr<Node> doubleWithLeftChild(std::shared_ptr<Node> pThis)
        {
            pRight = pRight->rotateWithLeftChild(pRight);
            return rotateWithRightChild(pThis);
        }//function

        virtual std::shared_ptr<Node> deleteFirst(std::shared_ptr<Node> pThis)
        {
            //check for leaf in pLeft
            if(pLeft->getHeight() != 0)
            {
                pLeft = pLeft->deleteFirst(pLeft);
                reCalcHeight();
                return pThis;
            }//if
            else
            {
                if(pNext.lock() != nullptr)
                    pNext.lock()->pPrev = std::shared_ptr<Branch>(nullptr);
                if(pRight->getHeight() != 0)
                    return pRight;
                else
                    return std::shared_ptr<Node>(new Leaf());
            }//else
            reCalcHeight();
        }//function

        virtual const std::shared_ptr<Node> getLeft() const
        {
            return pLeft;
        }//function

        virtual const std::shared_ptr<Node> getRight() const
        {
            return pRight;
        }//function

        virtual const std::shared_ptr<Node> getLeftOrLeaf(std::shared_ptr<Node> p) const
        {
            return pLeft;
        }//function

        virtual const std::shared_ptr<Node> getRightOrLeaf(std::shared_ptr<Node> p) const
        {
            return pRight;
        }//function

        void setNext(std::shared_ptr<Branch> pNew)
        {
            pNext = pNew;
        }//function

        void setPrev(std::shared_ptr<Branch> pNew)
        {
            pPrev = pNew;
        }//function

        std::shared_ptr<Branch> getNext() const
        {
            return pNext.lock();
        }//function

        std::shared_ptr<Branch> getPrev() const
        {
            return pPrev.lock();
        }//function
    };//class

    /**
     * @brief Self balancing binary search tree branch.
     */
    class Branch : public Node
    {
    private:
        T data;
        unsigned int height;

        //overload
        void setHeight(unsigned int h)
        {
            height = h;
        }//function

        virtual bool operator<(T& other) const
        {
            return data < other;
        }//operator

    public:
        Branch(T& data, std::shared_ptr<Branch> pPrev, std::shared_ptr<Branch> pNext,       std::shared_ptr<Node> pLeaf)
                :
            Node(pPrev, pNext, pLeaf, pLeaf),
            data(data),
            height(0)
        {}//constructor

        //overload
        const T& get() const
        {
            return data;
        }//function

        //overload
        unsigned int getHeight() const
        {
            return height;
        }//function

    };//class

    /**
     * @brief Self balancing binary search tree leaf.
     */
    class Leaf : public Node
    {
    private:
        //overload
        void setHeight(unsigned int h)
        {}//function

        bool operator<(T& other) const
        {
            throw NullPointerException("trying to compare with data of tree leaf");
        }//operator

    public:
        Leaf()
                :
            Node(
                    std::shared_ptr<Branch>(nullptr),
                    std::shared_ptr<Branch>(nullptr),
                    std::shared_ptr<Node>(nullptr),
                    std::shared_ptr<Node>(nullptr)
                )
        {}//constructor

        const T& get() const
        {
            throw NullPointerException("trying to access data of tree leaf");
        }//function

        //overload
        unsigned int getHeight() const
        {
            return 0;
        }//function

        std::shared_ptr<Node> rotateWithRightChild(std::shared_ptr<Node> pThis)
        {
            return pThis;
        }//function

        std::shared_ptr<Node> rotateWithLeftChild(std::shared_ptr<Node> pThis)
        {
            return pThis;
        }//function

        std::shared_ptr<Node> doubleWithRightChild(std::shared_ptr<Node> pThis)
        {
            return pThis;
        }//function

        std::shared_ptr<Node> doubleWithLeftChild(std::shared_ptr<Node> pThis)
        {
            return pThis;
        }//function
        
        std::shared_ptr<Node> insert(
                T& data,
                std::shared_ptr<Node> pThis,
                std::shared_ptr<Branch> pLastPrev,
                std::shared_ptr<Branch> pLastNext,
                std::shared_ptr<Branch>& pNew
            )
        {
            pNew = std::shared_ptr<Branch>(
                    new Branch(data, pLastPrev, pLastNext, pThis)
                );
            if(pLastPrev != nullptr)
                pLastPrev->setNext(pNew);
            if(pLastNext != nullptr)
                pLastNext->setPrev(pNew);
            assert(pNew->getRight() != nullptr);
            assert(pNew->getLeft() != nullptr);
            return pNew;
        }//function
        
        std::shared_ptr<Node> insertOrGet(
                T& data,
                std::shared_ptr<Node> pThis,
                std::shared_ptr<Branch> pLastPrev,
                std::shared_ptr<Branch> pLastNext,
                std::shared_ptr<Branch>& pNew
            )
        {
            return insert(data, pThis, pLastPrev, pLastNext, pNew);
        }//function

        const std::shared_ptr<Node> getLeft() const
        {
            throw NullPointerException("trying to access successor of tree leaf");
        }//function

        const std::shared_ptr<Node> getRight() const
        {
            throw NullPointerException("trying to access successor of tree leaf");
        }//function

        const std::shared_ptr<Node> getLeftOrLeaf(std::shared_ptr<Node> p) const
        {
            return p;
        }//function

        const std::shared_ptr<Node> getRightOrLeaf(std::shared_ptr<Node> p) const
        {
            return p;
        }//function
    };//class

    std::shared_ptr<Node> pRoot;
public:
    SearchTree()
            :
        pRoot(new Leaf())
    {}//constructor

    /** @brief Function to check if tree is empty */
    bool isEmpty() const
    {
        return pRoot->getHeight() == 0;
    }

    class Iterator
    {
    private:
        std::shared_ptr<Branch> pCurr;
    public:
        
        Iterator(std::shared_ptr<Branch> pCurr)
                :
            pCurr(pCurr)
        {}//constructor
        Iterator(const Iterator& xOther)
                :
            pCurr(xOther.pCurr)
        {}//constructor

        /**
         * @brief Get the content of the iterator.
         */
        const T& operator*() const 
        {
            if(pCurr == nullptr) 
                throw NullPointerException("trying to access iterator that has no element"); 
            return pCurr->get();
        }//function

        /**
         * @brief Get the content of the iterator.
         */
        const bool exists() const 
        {
            return pCurr != nullptr;
        }//function

        /**
         * @brief Get the content of the iterator.
         */
        const T& operator->() const
        {
            if(pCurr == nullptr) 
                throw NullPointerException("trying to access iterator that has no element"); 
            return pCurr->get();
        }//function

        /**
         * @brief Move to the next list element.
         */
        void operator++()
        {
            if(pCurr != nullptr) 
                pCurr = pCurr->next();
        }//function
    };//class

    /**
     * @brief Function to insert data.
     * @returns the new tree element.
     */
    Iterator insert(T& data)
    {
        std::shared_ptr<Branch> pRet = nullptr;
        pRoot = pRoot->insert(data, pRoot, nullptr, nullptr, pRet);
        assert(pRoot != nullptr);
        assert(pRet != nullptr);
        return Iterator(pRet);
    }

    /**
     * @brief Function to insert a new element if not already present.
     * @returns the element with the given value if present; the new element otherwise.
     */
    Iterator insertOrGet(T& data)
    {
        std::shared_ptr<Branch> pRet = nullptr;
        pRoot = pRoot->insertOrGet(data, pRoot, nullptr, nullptr, pRet);
        assert(pRoot != nullptr);
        assert(pRet != nullptr);
        return Iterator(pRet);
    }

    /**
     * @brief Delete the leftmost element of the tree.
     */
    void deleteFirst()
    {
        pRoot = pRoot->deleteFirst(pRoot);
    }//function

    /** @brief returns the leftmost element of the Tree. */
    const T& first() const
    {
        assert(!isEmpty());
        std::shared_ptr<Node> pCurr = pRoot;
        while(pCurr->getHeight() > 1)
        {
            if(pCurr->getLeft()->getHeight() > 0)
                pCurr = pCurr->getLeft();
            else
                pCurr = pCurr->getRight();
        }
        return pCurr->get();
    }//function

    /** @brief returns the root element of the Tree. */
    const T& root() const
    {
        assert(!isEmpty());
        return pRoot->get();
    }//function
};//class

#endif