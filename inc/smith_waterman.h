#ifndef SMITH_WATERMAN
#define SMITH_WATERMAN

#include "cppModule.h"
#include "alignment.h"

void exportSMW();

class SMW: public CppModule
{
    //overload
    std::shared_ptr<Container> execute(ContainerVector vpInput);

    //overload
    ContainerVector getInputType() const;

    //overload
    std::shared_ptr<Container> getOutputType() const;

    std::string getName() const
    {
        return "SMW";
    }
};//class

#endif