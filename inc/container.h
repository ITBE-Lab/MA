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


class Container
{
public:
    virtual ContainerType getType(){return ContainerType::unknown;}
    
    virtual std::string getTypeInfo();
};//class

/* function to export this module to boost python */
void exportContainer();
#endif