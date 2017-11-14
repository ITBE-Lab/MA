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

/** 
 * @brief Function to export Container to boost python.
 * @ingroup export
 */
void exportContainer();
#endif