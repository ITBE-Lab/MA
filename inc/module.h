/** 
 * @file module.h
 * @brief Implements a abstract Module class.
 * @author Markus Schmidt
 */
#ifndef MODULE_H
#define MODULE_H

#include <container.h>
#include <memory>
#include <Python.h>
#include "exception.h"
#include <iostream>


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
    std::shared_ptr<Container> saveExecute(std::vector<std::shared_ptr<Container>> pInput)
    {
        std::cout << "TODO: save exec disabled" << std::endl;
        return execute(pInput);
    }//function

    FutureContainer promiseMe(std::vector<FutureContainer> vInput)
    {
        //TODO: check types
        return FutureContainer(this, getOutputType(), vInput);
    }//function
};

/**
 * @brief Exposes the Module class to boost python.
 */
void exportModule();

#endif