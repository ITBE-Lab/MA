/** 
 * @file module.h
 * @brief Implements a abstract Module class.
 * @author Markus Schmidt
 */
#ifndef MODULE_H
#define MODULE_H

#include "container.h"
#include <memory>
#include <Python.h>
#include <iostream>

class Pledge;


bool typeCheck(
        std::shared_ptr<Container> pData, 
        ContainerType xExpected
    );

bool typeCheck(
        std::vector<std::shared_ptr<Container>> vData, 
        std::vector<ContainerType> vExpected
    );

/**
 * @brief Abstract class intended for the implementaiton of various algorithms.
 * @details
 * All computing on data should inherit from this class
 */
class Module
{
public:
    /**
     * @brief Execute the implemented algorithm
     * @details
     * Expects the given containers to have the correct types.
     */
    virtual std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> pInput)
    {
        return nullptr;
    }

    /**
     * @brief The expected input types.
     * @details
     * Used for type checking the inputs before calling execute.
     */
    virtual std::vector<ContainerType> getInputType()
    {
        return std::vector<ContainerType>{ContainerType::nothing};
    }

    /**
     * @brief The expected output type.
     * @details
     * Used for type checking weather the module returns expected data.
     */
    virtual ContainerType getOutputType()
    {
        return ContainerType::nothing;
    }

    /**
     * @brief Execute the implemented algorithm
     * @details
     * Internally calls execute after checking the input types.
     * Also checks the result returned by Execute.
     */
    std::shared_ptr<Container> saveExecute(std::vector<std::shared_ptr<Container>> vInput)
    {
        if(!typeCheck(vInput, getInputType()))
            throw new ModuleIO_Exception("Input type and expected input type did not match.");
        std::shared_ptr<Container> pRet = execute(vInput);
        if(!typeCheck(pRet, getOutputType()))
            throw new ModuleIO_Exception("Module produced output of wrong type.");
        return pRet;
    }//function

    Pledge promiseMe(std::vector<std::shared_ptr<Pledge>> vInput);
};

/**
 * @brief Abstract class intended to hold promises to data objects used by Modules.
 */
class Pledge : public Container
{
private:
    Module* pledger;
    std::shared_ptr<Container> content;
    ContainerType type;
    unsigned int iLastCallNum = 0;
    std::vector<std::shared_ptr<Pledge>> vPredecessors;
public:
    Pledge(
            ContainerType type
        )
            :
        pledger(),
        content(),
        type(type),
        vPredecessors()
    {}//constructor

    Pledge(
            Module* pledger,
            ContainerType type,
            std::vector<std::shared_ptr<Pledge>> vPredecessors
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
        if(iLastCallNum < iCallNum)
        {
            if(pledger == nullptr)
                throw ModuleIO_Exception("No pledger known");
            iLastCallNum = iCallNum;
            std::vector<std::shared_ptr<Container>> vInput;
            for(std::shared_ptr<Pledge> pFuture : vPredecessors)
                vInput.push_back(pFuture->fullfill(iCallNum));
            content = pledger->execute(vInput);
            assert(content->getType());
        }//if
        return content;
    }//function
    
    std::shared_ptr<Container> get()
    {
        return fullfill(iLastCallNum);
    }//function

    //overload
    ContainerType getType()
    {
        return type;
    }//function
};//class

/**
 * @brief Exposes the Module class to boost python.
 */
void exportModule();

#endif