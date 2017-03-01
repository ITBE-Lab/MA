#ifndef ALIGNER_H
#define ALIGNER_H

#include "modules/module.h"
#include "data/container.h"
#include <list>
#include <memory>
#include <boost/python.hpp>
#include <Python.h>
#include <iostream>

class Aligner
{
private:
    std::list<Module> lModules;
    std::shared_ptr<Container> pCurrent;

public:
    Aligner()
        :
        lModules(), pCurrent()
    {}

    void setData(std::shared_ptr<Container> pC){pCurrent = pC;}
    void addModule(Module xM){lModules.push_back(xM);std::cout << "added Module" << std::endl;}
    void step(){if(lModules.empty())return; pCurrent = lModules.front().execute(pCurrent); lModules.pop_front();}
    void steps(){while(!lModules.empty())step();}

    void test()
    {
        addModule(Printer());
        step();
    }

};

#endif