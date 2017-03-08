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
    virtual std::shared_ptr<Container> execute(std::shared_ptr<Container> pInput){return pInput;}

    virtual std::shared_ptr<Container> getInputType(){return std::shared_ptr<Container>(new DummyContainer(ContainerType::nothing));}

    virtual std::shared_ptr<Container> getOutputType(){return std::shared_ptr<Container>(new DummyContainer(ContainerType::nothing));}

    std::shared_ptr<Container> saveExecute(std::shared_ptr<Container> pInput)
    {
        if(getInputType() == nullptr) 
                throw ModuleIO_Exception("expected input of module is a nullptr - use DummyContainer(nothing) instead");
        if(pInput == nullptr)
                throw ModuleIO_Exception("given input is a nullptr");
        if(getInputType()->sameTypeAs(pInput))
        {
            auto pOutput = execute(pInput);
            if(getOutputType() == nullptr) 
                throw ModuleIO_Exception("expected output of module is a nullptr - use DummyContainer(nothing) instead");
            if(pOutput == nullptr)
                throw ModuleIO_Exception("module returned a nullptr");
            if(getOutputType()->sameTypeAs(pOutput))
                return pOutput;
            else
                throw ModuleIO_Exception(
                        std::string("wrong output type. Expected: ").append(getOutputType()->getTypeInfo())
                        .append(", but got: ").append(pOutput->getTypeInfo())
                        .c_str() );
        }//if
        else
            throw ModuleIO_Exception(
                        std::string("wrong input type. Expected: ").append(getInputType()->getTypeInfo())
                        .append(", but got: ").append(pInput->getTypeInfo())
                        .c_str() );
    }//function
};

class Printer : public Module
{
public:
    Printer(){}
    // overwrite the execute frunction from Module
    std::shared_ptr<Container> execute(std::shared_ptr<Container> pInput) override
    {
        if(pInput != nullptr)
            pInput->print();
        else
            std::cout << "nullptr" << std::endl;
        return pInput;
    }
    std::shared_ptr<Container> getInputType()
    {
        return std::shared_ptr<Container>(new DummyContainer(ContainerType::any));
    }
    std::shared_ptr<Container> getOutputType()
    {
        return std::shared_ptr<Container>(new DummyContainer(ContainerType::any));
    }
};

//function called in order to export this module
void exportModule();

#endif