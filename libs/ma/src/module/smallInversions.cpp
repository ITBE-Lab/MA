/**
 * @file smallInversions.cpp
 * @author Markus Schmidt
 */
#include "ma/module/smallInversions.h"


using namespace libMA;
using namespace libMS;


#ifdef WITH_PYTHON

void exportSmallInversions( py::module& rxPyModuleId )
{
    // export the SmallInversions class
    exportModule<SmallInversions>( rxPyModuleId, "SmallInversions" );
} // function
#endif