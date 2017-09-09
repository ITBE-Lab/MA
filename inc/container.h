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
    *   return the type of the container as enum
    *   used by Module for type checking the inputs
    */
    virtual ContainerType getType(){return ContainerType::unknown;}
    
    virtual std::string getTypeInfo();
};//class

/** function to export this module to boost python */
void exportContainer();
#endif