#ifndef ALIGNER_H
#define ALIGNER_H

#include "module.h"
#include "container.h"
#include "segmentation.h"
#include "FM_index.h"
#include "container.h"
#include "sequence.h"
#include <list>
#include <memory>
#include <boost/python.hpp>
#include <Python.h>
#include <iostream>

/*
*	implements the alignment
*/


class Aligner
{
private:
    std::list<std::shared_ptr<Module>> lpModules;
    std::shared_ptr<Container> pCurrent;

public:
    Aligner()
        :
        lpModules(), pCurrent(new DummyContainer(ContainerType::nothing))
    {}

    void setData(std::shared_ptr<Container> pC){pCurrent = pC;}
   
    void addModule(Module xM)
    {
        lpModules.push_back(std::make_shared<Module>(xM));
    }
    bool stepPossible() 
    {
        if(lpModules.empty())
            return false; 
        if(pCurrent == nullptr)
            return false;
        if(lpModules.front() == nullptr)
            return false;
        return true; 
    }
    void step()
    {
        if(!stepPossible())
            return; 
        pCurrent = lpModules.front()->saveExecute(pCurrent); 
        lpModules.pop_front();
    }
    void steps()
    {
        while(stepPossible())
            step();
    }

    void test()
    {
        addModule(Printer());
        steps();
    }

};

#endif