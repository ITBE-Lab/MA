/**
 * @file smallInversions.cpp
 * @author Markus Schmidt
 */
#include "ma/module/smallInversions.h"


using namespace libMA;
using namespace libMS;


#ifdef WITH_PYTHON

void exportSmallInversions( libMS::SubmoduleOrganizer& xOrganizer )
{
    // export the SmallInversions class
    exportModule<SmallInversions>( xOrganizer, "SmallInversions" );
} // function
#endif