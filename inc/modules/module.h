#ifndef MODULE_H
#define MODULE_H

#include <data/container.h>
#include <memory>
#include <boost/python.hpp>
#include <Python.h>
#include <iostream>

class Module
{
public:
    virtual std::shared_ptr<Container> execute(std::shared_ptr<Container> pInput);

};

class Printer : public Module
{
     std::shared_ptr<Container> execute(std::shared_ptr<Container> pInput) override;
};

#endif