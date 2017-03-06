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
        std::cout << "Module.saveExecute()" << std::endl;
        if(this == nullptr) 
                throw ModuleIO_Exception("this is a nullptr WTF !?!");
        if(getInputType() == nullptr) 
                throw ModuleIO_Exception("expected input of module is a nullptr - use DummyContainer(nothing) instead");
        std::cout << "Module.saveExecute() - 0" << std::endl;
        if(pInput == nullptr)
                throw ModuleIO_Exception("given input is a nullptr");
        std::cout << "Module.saveExecute() - 1" << std::endl;
        if(getInputType()->sameTypeAs(pInput))
        {
        std::cout << "Module.saveExecute() - 2" << std::endl;
            auto pOutput = execute(pInput);
            if(getOutputType() == nullptr) 
                throw ModuleIO_Exception("expected output of module is a nullptr - use DummyContainer(nothing) instead");
            if(pOutput == nullptr)
                throw ModuleIO_Exception("module returned a nullptr");
            if(getOutputType()->sameTypeAs(pOutput))
                return pOutput;
            else
                throw ModuleIO_Exception("wrong output type");
        }
        else
            throw ModuleIO_Exception("wrong input type");
    }
};

class Printer : public Module
{
public:
    Printer(){std::cout << "Printer()" << std::endl;}
    // printing deconstructor to verify the livetime of the object
    ~Printer(){std::cout << "~Printer()" << std::endl;}
    // overwrite the execute frunction from Module
    std::shared_ptr<Container> execute(std::shared_ptr<Container> pInput) override
    {
        std::cout << "Printer.execute()" << std::endl; 
        return pInput;
    }
    std::shared_ptr<Container> getInputType()
    {
        std::cout << "Printer.getInputType()" << std::endl;
        return std::shared_ptr<Container>(new DummyContainer(ContainerType::nothing));
    }
    std::shared_ptr<Container> getOutputType()
    {
        return std::shared_ptr<Container>(new DummyContainer(ContainerType::nothing));
    }
};

//function called in order to export this module
void exportModule();

#endif