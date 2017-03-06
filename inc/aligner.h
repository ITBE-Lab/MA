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
    /* pyhton will take care of the deletion the the given pointer */
    void addModule(std::shared_ptr<Module> pM)
    {
        lpModules.push_back(pM);
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
        std::cout << "Aligner.step()" << std::endl;
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
        addModule(std::shared_ptr<Module>(new Printer()));
        steps();
    }

};

#endif