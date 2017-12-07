/**
 * @file container.h
 * @brief Implements Container.
 * @author Markus Schmidt
 */

#ifndef CONTAINER_H
#define CONTAINER_H

#include <vector>
#include "exception.h"
#include "iterableConverter.h"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

/**
 * @defgroup container
 * @brief All classes containing data elements.
 * @details
 * All classes that store data should inherit from Container.
 * These classes are then added to this Group.
 */

/**
 * @brief Abstract class intended to hold data objects used by @ref CppModule "modules".
 * @details
 * All classes containing data should inherit from this class.
 * @ingroup container
 */
class Container
{
public:
    /** 
    * @returns the type of the container as int.
    * @brief Used by @ref CppModule "module" for type checking its Inputs.
    */
    volatile bool canCast(std::shared_ptr<Container> c) const
    {
        return true;
    }//function

    virtual std::string getTypeName() const
    {
        return "Container";
    }//function

    virtual std::shared_ptr<Container> getType() const
    {
        return std::shared_ptr<Container>(new Container());
    }//function

    /**
     * @details
     * used by ExecOnVec to sort containers in order
     */
    virtual bool smaller(const std::shared_ptr<Container> other) const
	{
		return true;
	}// operator
};//class


class Nil : public Container
{
public:
    /** 
    * @returns the type of the container as int.
    * @brief Used by @ref CppModule "module" for type checking its Inputs.
    */
    bool canCast(std::shared_ptr<Container> c) const
    {
        return false;
    }//function

    std::string getTypeName() const
    {
        return "Nil";
    }//function

    std::shared_ptr<Container> getType() const
    {
        return std::shared_ptr<Container>(new Nil());
    }//function
};//class

class ContainerVector:
    public Container,
    public std::vector<std::shared_ptr<Container>>
{
private:
    std::shared_ptr<Container> contentType;
public:
    using vector::vector;
    ContainerVector()
        :
        contentType(new Container())
    {}//container vector

    ContainerVector(std::shared_ptr<Container> contentType)
        :
        contentType(contentType)
    {}//container vector

    ContainerVector(const std::shared_ptr<ContainerVector> pOther)
            :
        vector(*pOther),
        contentType(pOther->contentType)
    {}//container vector

    ContainerVector(std::shared_ptr<std::vector<std::shared_ptr<Container>>> pContent)
            :
        vector(*pContent),
        contentType(pContent->front()->getType())
    {}//container vector

    ContainerVector(std::shared_ptr<Container> contentType, size_t numElements)
        :
        vector(numElements),
        contentType(contentType)
    {}//container vector

    //overload
    bool canCast(std::shared_ptr<Container> c) const
    {
        std::shared_ptr<ContainerVector> casted = std::dynamic_pointer_cast<ContainerVector>(c);
        if(casted == nullptr)
            return false;
        return casted->contentType->getType()->canCast(contentType->getType());
    }//function

    //overload
    std::string getTypeName() const
    {
        return "ContainerVector";
    }//function

    //overload
    std::shared_ptr<Container> getType() const
    {
        return std::shared_ptr<Container>(new ContainerVector(contentType));
    }//function

    std::vector<std::shared_ptr<Container>> get()
    {
        return *this;
    }//function

    void push_back_boost(std::shared_ptr<Container> c)
    {
        push_back(c);
    }//function

};//class

/** 
 * @brief Function to export Container to boost python.
 * @ingroup export
 */
void exportContainer();
#endif