#ifndef SMITH_WATERMAN
#define SMITH_WATERMAN

#include "cppModule.h"
#include "alignment.h"

void exportSMW();

class SMW: public CppModule
{
    //overload
    std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> vpInput);

    //overload
    std::vector<ContainerType> getInputType();

    //overload
    ContainerType getOutputType();
};//class

#endif