/**
 * @file smallInversions.cpp
 * @author Markus Schmidt
 */


#include "module/smallInversions.h"


using namespace libMA;


#ifdef WITH_PYTHON

void exportSmallInversions( py::module& rxPyModuleId )
{
    // export the SmallInversions class
    exportModule<SmallInversions>( rxPyModuleId, "SmallInversions" );
} // function
#endif