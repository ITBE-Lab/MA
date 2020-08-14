/**
 * @file count_k_mers.cpp
 * @author Markus Schmidt
 */
#include "msv/module/count_k_mers.h"
#include <pybind11/stl.h>

using namespace libMSV;
using namespace libMS;

#ifdef WITH_PYTHON


void exportCountKMers( libMS::SubmoduleOrganizer& xOrganizer )
{
    py::class_<KMerCounter, Container, std::shared_ptr<KMerCounter>>( xOrganizer.container( ), "KMerCounter" )
        .def( py::init<nucSeqIndex>( ) );
    exportModule<GetKMerCounter, nucSeqIndex>( xOrganizer, "GetKMerCounter" );
    exportModule<KMerCounterModule>( xOrganizer, "KMerCounterModule" );
    exportModule<KMerCountFilterModule, nucSeqIndex>( xOrganizer, "KMerCountFilterModule" );

    py::class_<HashCounter, Container, std::shared_ptr<HashCounter>>( xOrganizer.container( ), "HashCounter" )
        .def( py::init<>( ) );
    exportModule<MMCounterModule>( xOrganizer, "MMCounterModule" );
    exportModule<MMFilteredSeeding, nucSeqIndex>( xOrganizer, "MMFilteredSeeding" );
} // function
#endif