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
    virtual std::shared_ptr<Container> getInputType(){return std::shared_ptr<Container>(new Container());}
    virtual std::shared_ptr<Container> getOutputType(){return std::shared_ptr<Container>(new Container());}
    std::shared_ptr<Container> saveExecute(std::shared_ptr<Container> pInput)
    {
        if(getInputType()->sameTypeAs(pInput))
        {
            auto pOutput = execute(pInput);
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
};

//function called in order to export this module
void exportModule();

#endif