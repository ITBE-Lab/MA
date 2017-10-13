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

/**
 * @brief Used to describe the type of the Container.
 * @details
 * required by CppModule for type checking the inputs before casting and using them.
 */
enum ContainerType{
    fM_index,
    nucSeq,
    alignment,
    packedNucSeq,
    segmentList,
    segment,
    seed,
    seeds,
    seedsVector,
    sa_interval,
    unknown,
    nothing,
    any,
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