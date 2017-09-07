#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <boost/python.hpp>
#include "container.h"

class Alignment : public Container
{
public:
    ContainerType getType(){return ContainerType::alignment;}
};//class

/* function to export this module to boost python */
void exportAlignment();
#endif