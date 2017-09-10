#ifndef CONTAINER_H
#define CONTAINER_H

#include <vector>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <iostream>

enum ContainerType{
    fM_index,
    nucSeq,
    alignment,
    packedNucSeq,
    segmentList,
    segment,
    stripOfConsideration,
    stripOfConsiderationList,
    vector,
    sa_interval,
    unknown,
    nothing,
    any,
};//enum

/**
*   all classes containing data should inherit from this class
*/
class Container
{
public:
    /** 
    *   returns the type of the container as enum.
    *   Used by Module for type checking its Inputs.
    */
    virtual ContainerType getType(){return ContainerType::unknown;}
};//class

/** function to export this module to boost python */
void exportContainer();
#endif