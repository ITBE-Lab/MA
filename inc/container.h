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
 * @brief Used to describe the type of the Container.
 * @details
 * required by CppModule for type checking the inputs before casting and using them.
 */
enum ContainerType{
    fM_index,//0
    nucSeq,//1
    alignment,//2
    packedNucSeq,//3
    segmentList,//4
    segment,//5
    seed,//6
    seeds,//7
    seedsVector,//8
    sa_interval,//9
    unknown,//10
    nothing,//11
    any,//12
};//enum

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
    * @returns the type of the container as enum.
    * @brief Used by @ref CppModule "module" for type checking its Inputs.
    */
    virtual ContainerType getType(){return ContainerType::unknown;}
};//class

/** 
 * @brief Function to export Container to boost python.
 * @ingroup export
 */
void exportContainer();
#endif