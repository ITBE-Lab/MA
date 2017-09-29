/** 
 * @file container.h
 * @brief Implements Container.
 * @author Markus Schmidt
 */

#ifndef CONTAINER_H
#define CONTAINER_H

#include <cmath>
#include <boost/python.hpp>
#include <vector>

/**
 * @brief Used to describe the type of the Container.
 * @details
 * required by module for type checking the inputs before casting and using them.
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
 * @brief Abstract class intended to hold data objects used by Modules.
 * @details
 * All classes containing data should inherit from this class.
 */
class Container
{
public:
    /** 
    * @returns the type of the container as enum.
    * @brief Used by Module for type checking its Inputs.
    */
    virtual ContainerType getType(){return ContainerType::unknown;}
};//class

/** 
 * @brief Function to export Container to boost python.
 */
void exportContainer();
#endif