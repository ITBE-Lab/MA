/**
 * @file sweepSvJumps.cpp
 * @author Markus Schmidt
 */
#include "msv/module/sweepSvJumps.h"
#include "ms/util/pybind11.h"

using namespace libMSV;
using namespace libMS;


#ifdef WITH_PYTHON

#include "ms/container/sv_db/py_db_conf.h"

bool AbstractFilter::bSilent = false;

void exportSweepSvJump( libMS::SubmoduleOrganizer& xOrganizer )
{
    py::class_<GenomeSection, Container, std::shared_ptr<GenomeSection>>( xOrganizer.container( ), "GenomeSection" )
        .def_readwrite( "start", &GenomeSection::iStart )
        .def_readwrite( "size", &GenomeSection::iSize );

    py::bind_vector<std::vector<std::shared_ptr<SvCall>>>( xOrganizer.util( ), "SvCallVec", "docstr" );

    py::class_<CompleteBipartiteSubgraphClusterVector,
               Container,
               std::shared_ptr<CompleteBipartiteSubgraphClusterVector>>( xOrganizer.container( ),
                                                                         "CompleteBipartiteSubgraphClusterVector" )
        .def( py::init<>( ) )
        .def_readwrite( "content", &CompleteBipartiteSubgraphClusterVector::vContent );

    exportModule<GenomeSectionFactory, std::shared_ptr<Pack>>( xOrganizer, "GenomeSectionFactory" );


    exportModule<CompleteBipartiteSubgraphSweep<DBCon>, int64_t>(
        xOrganizer, "CompleteBipartiteSubgraphSweep", []( auto&& x ) {
            x.def_readonly( "time_init", &CompleteBipartiteSubgraphSweep<DBCon>::dInit )
                .def_readonly( "time_complete_while", &CompleteBipartiteSubgraphSweep<DBCon>::dOuterWhile )
                .def_readonly( "time_inner_while", &CompleteBipartiteSubgraphSweep<DBCon>::dInnerWhile );
        } );

    exportModule<ExactCompleteBipartiteSubgraphSweep>( xOrganizer, "ExactCompleteBipartiteSubgraphSweep" );
    exportModule<FilterLowSupportShortCalls>( xOrganizer, "FilterLowSupportShortCalls" );
    exportModule<FilterLowScoreCalls>( xOrganizer, "FilterLowScoreCalls" );
    exportModule<FilterDiagonalLineCalls>( xOrganizer, "FilterDiagonalLineCalls" );
    exportModule<FilterFuzzyCalls>( xOrganizer, "FilterFuzzyCalls" );
    exportModule<ComputeCallAmbiguity>( xOrganizer, "ComputeCallAmbiguity" );

    py::class_<AbstractFilter>( xOrganizer.util( ), "AbstractFilter" )
        .def_readwrite_static( "silent", &AbstractFilter::bSilent );
} // function
#endif