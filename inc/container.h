#ifndef CONTAINER_H
#define CONTAINER_H

#include <vector>
#include <boost/python.hpp>
#include <iostream>

enum ContainerType{
    fM_index,
    nucSeq,
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