/** 
 * @file module.h
 * @brief Implements a abstract Module class.
 * @author Markus Schmidt
 */
#ifndef MODULE_H
#define MODULE_H

#include "container.h"
#include "debug.h"
#include <memory>
#include <Python.h>
#include <iostream>
#include <boost/python/list.hpp>

class Pledge;

/**
 * @brief Checks weather the given types are equal.
 * @details
 * This requires a function since we have types like any of none.
 */
bool typeCheck(
        ContainerType xData, 
        ContainerType xExpected
    );

/**
 * @brief Checks weather the given container has the expected type.
 * @details
 * This requires a function since we have types like any of none.
 */
bool typeCheck(
        std::shared_ptr<Container> pData, 
        ContainerType xExpected
    );

/**
 * @brief Checks weather given containers have the expected types.
 * @details
 * This requires a function since we have types like any of none.
 */
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

    /**
     * @brief Make this module promise to execute it's function on the provided data.
     * @details
     * see Pledge for more information about the computational graph that can be build using 
     * promiseMe and Pledge.
     */
    std::shared_ptr<Pledge> promiseMe(std::vector<std::shared_ptr<Pledge>> vInput);
};

/**
 * @brief Abstract class intended to hold promises to data objects used by Modules.
 * @details
 * Us this class to set up a computational graph.
 * The graphs nodes are Modules (in python or cpp) that perform computational tasks.
 * Each node takes a vector of containers as input.
 * This container vector was created using the vector of Pledges that was given to the Module.
 * Then the Module returns a Container that fits the pledge returned by the Modules 
 * promiseMe function.
 *
 * @note The advantage of the computational graph is that there is no unnecessary 
 * jumping between python and cpp code.
 */
class Pledge : public Container
{
private:
    Module* pledger;
    boost::python::object py_pledger;
    std::shared_ptr<Container> content;
    ContainerType type;
    unsigned int iLastCallNum = 0;
    static unsigned int iCurrentCallNum;
    std::vector<std::shared_ptr<Pledge>> vPredecessors;
public:
    /**
     * @brief Create a new pledge without a module giving the pledge.
     * @details
     * This means that this Pledge has to be fullfilled by calling set manually.
     */
    Pledge(
            ContainerType type
        )
            :
        pledger(nullptr),
        py_pledger(),
        content(),
        type(type),
        vPredecessors()
    {}//constructor

    /**
     * @brief Create a new pledge with a cpp module responsible to fullfill it.
     * @details
     * This means that this Pledge can be automatically fullfilled by the given module.
     */
    Pledge(
            Module* pledger,
            ContainerType type,
            std::vector<std::shared_ptr<Pledge>> vPredecessors
        )
            :
        pledger(pledger),
        py_pledger(),
        content(),
        type(type),
        vPredecessors(vPredecessors)
    {}//constructor

    /**
     * @brief Create a new pledge with a Python module responsible to fullfill it.
     * @details
     * This means that this Pledge can be automatically fullfilled by the given module.
     */
    Pledge(
            boost::python::object py_pledger,
            ContainerType type,
            std::vector<std::shared_ptr<Pledge>> vPredecessors
        )
            :
        pledger(),
        py_pledger(py_pledger),
        content(),
        type(type),
        vPredecessors(vPredecessors)
    {}//constructor

    /**
     * @brief Manually fullfill this pledge.
     * @note It would be better to have a "constant" or "variable" class instead of this function.
     */
    void set(std::shared_ptr<Container> c)
    {
        iLastCallNum = iCurrentCallNum+1;
        content = c;
    }//function

    /**
     * @brief automatically fullfill this pledge.
     * @details
     * Checks weather the pledge has already been fullfilled for this call.
     * If necessary, uses the python or cpp module to fullfill the pledge,
     * and returns the respective container.
     */
    std::shared_ptr<Container> fullfill(unsigned int iCallNum)
    {
        if(iLastCallNum < iCallNum)
        {
            if(pledger == nullptr && py_pledger.is_none())
                throw ModuleIO_Exception("No pledger known");
            iLastCallNum = iCallNum;
            if(pledger != nullptr)
            {
                std::vector<std::shared_ptr<Container>> vInput;
                for(std::shared_ptr<Pledge> pFuture : vPredecessors)
                    vInput.push_back(pFuture->fullfill(iCallNum));
                content = (std::shared_ptr<Container>)pledger->execute(vInput);
            }//if
            else
            {
                boost::python::list vInput;
                for(std::shared_ptr<Pledge> pFuture : vPredecessors)
                    vInput.append(pFuture->fullfill(iCallNum));
                /*
                 * here we jump to python code to call a function and resume the cpp code 
                 * once python is done...
                 */ 
                content = boost::python::extract<std::shared_ptr<Container>>(py_pledger.attr("execute")(vInput));
            }//else
            assert(typeCheck(content, type));
        }//if
        return content;
    }//function

    /**
     * @brief Get the promised container.
     * @details
     * Tries to automatically fullfill the pledge, if necessary.
     */
    std::shared_ptr<Container> get()
    {
        return fullfill(iCurrentCallNum);
    }//function
    
    /**
     * @brief Get the promised container.
     * @details
     * Invalidates the content of all pledges that have not been manually set.
     * Then fullfills all necessary pledges to return the promised content of this pledge.
     * If multiple alignments shall be done, this can be used to do so.
     */
    std::shared_ptr<Container> next()
    {
        return fullfill(++iCurrentCallNum);
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