#define LINESWEEP
#ifndef LINESWEEP
#define LINESWEEP

#include "module.h"

class LineSweep: public Module
{
public:

	LineSweep(){}//constructor

	std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> pInput);

    std::vector<ContainerType> getInputType();

    std::vector<ContainerType> getOutputType();
};//class

void exportLinesweep();

#endif