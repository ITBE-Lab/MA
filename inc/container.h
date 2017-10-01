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
#include "exception.h"

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
 * @brief Abstract class intended to hold promises to data objects used by Modules.
 */
class FutureContainer
{
private:
    Module* pledger;
    std::shared_ptr<Container> content;
    ContainerType type;
    unsigned int iLastCallNum = 0;
    std::vector<FutureContainer> vPredecessors;
public:
    FutureContainer(
            ContainerType type
        )
            :
        pledger(),
        content(),
        type(type),
        vPredecessors()
    {}//constructor

    FutureContainer(
            Module* pledger,
            ContainerType type,
            std::vector<FutureContainer> vPredecessors
        )
            :
        pledger(pledger),
        content(),
        type(type),
        vPredecessors(vPredecessors)
    {}//constructor

    void set(std::shared_ptr<Container> c, unsigned int iCallNum)
    {
        iLastCallNum = iCallNum;
        content = c;
    }//function

    std::shared_ptr<Container> fullfill(unsigned int iCallNum)
    {
        if(pledger == nullptr)
            throw ModuleIO_Exception("No pledger and no data for future provided");
        if(iLastCallNum < iCallNum)
        {
            iLastCallNum = iCallNum;
            std::vector<std::shared_ptr<Container>> vInput;
            for(FutureContainer xFuture : vPredecessors)
                vInput.push_back(xFuture.fullfill(iCallNum));
            content = pledger.execute(vInput);
        }//if
        if(content == nullptr)
            throw ModuleIO_Exception("pledger did not hold his promise");
        return content;
    }//function
    
    std::shared_ptr<Container> get()
    {
        return fullfill(iLastCallNum);
    }//function
};//class

/** 
 * @brief Function to export Container to boost python.
 */
void exportContainer();
#endif