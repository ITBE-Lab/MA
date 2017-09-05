#ifndef MODULE_H
#define MODULE_H

#include <container.h>
#include <memory>
#include <boost/python.hpp>
#include <Python.h>
#include "exception.h"
#include <iostream>

/*
    Module is an abstract class, but boost python requires the functions to exist,
    however we can remove the ability to initialize this class for python
*/
class Module
{
public:
    virtual std::vector<std::shared_ptr<Container>> execute(
            std::vector<std::shared_ptr<Container>> pInput
        )
    {
        return pInput;
    }

    virtual std::ContainerType<ContainerType> getInputType()
    {
        return std::ContainerType<ContainerType>(ContainerType::nothing);
    }

    virtual std::vector<ContainerType> getOutputType()
    {
        return std::vector<Container>(ContainerType::nothing);
    }

    std::vector<std::shared_ptr<Container>> saveExecute(
            std::vector<std::shared_ptr<Container>> pInput
        )
    {
        std::cout << "TODO: save exec disabled" << std::endl;
        return execute(pInput);
    }//function
};

class Printer : public Module
{
public:
    Printer(){}
    // overwrite the execute frunction from Module
    std::vector<std::shared_ptr<Container>> execute(
            std::vector<std::shared_ptr<Container>> pInput
        ) override
    {
        
        std::cout << "TODO: printer disabled" << std::endl;
        return pInput;
    }

    virtual std::ContainerType<ContainerType> getInputType()
    {
        return std::ContainerType<ContainerType>(ContainerType::any);
    }

    virtual std::vector<ContainerType> getOutputType()
    {
        return std::vector<Container>(ContainerType::any);
    }
};

//function called in order to export this module
void exportModule();

#endif