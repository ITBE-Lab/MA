/**
 * @file sweepSvJumps.cpp
 * @author Markus Schmidt
 */
#include "module/sweepSvJumps.h"
#include "util/pybind11.h"

using namespace libMSV;


#ifdef WITH_PYTHON

#include "container/sv_db/py_db_conf.h"

void exportSweepSvJump( py::module& rxPyModuleId )
{
    py::class_<GenomeSection, Container, std::shared_ptr<GenomeSection>>( rxPyModuleId, "GenomeSection" )
        .def_readwrite( "start", &GenomeSection::iStart )
        .def_readwrite( "size", &GenomeSection::iSize );

    py::bind_vector<std::vector<std::shared_ptr<SvCall>>>( rxPyModuleId, "SvCallVec", "docstr" );

    py::class_<CompleteBipartiteSubgraphClusterVector,
               Container,
               std::shared_ptr<CompleteBipartiteSubgraphClusterVector>>( rxPyModuleId,
                                                                         "CompleteBipartiteSubgraphClusterVector" )
        .def( py::init<>( ) )
        .def_readwrite( "content", &CompleteBipartiteSubgraphClusterVector::vContent );

    exportModule<GenomeSectionFactory, std::shared_ptr<Pack>>( rxPyModuleId, "GenomeSectionFactory" );


    exportModule<CompleteBipartiteSubgraphSweep<DBCon>, int64_t>(
        rxPyModuleId, "CompleteBipartiteSubgraphSweep", []( auto&& x ) {
            x.def_readonly( "time_init", &CompleteBipartiteSubgraphSweep<DBCon>::dInit )
                .def_readonly( "time_complete_while", &CompleteBipartiteSubgraphSweep<DBCon>::dOuterWhile )
                .def_readonly( "time_inner_while", &CompleteBipartiteSubgraphSweep<DBCon>::dInnerWhile );
        } );

    exportModule<ExactCompleteBipartiteSubgraphSweep<DBCon>>( rxPyModuleId, "ExactCompleteBipartiteSubgraphSweep" );
    exportModule<FilterLowSupportShortCalls>( rxPyModuleId, "FilterLowSupportShortCalls" );
    exportModule<FilterLowScoreCalls>( rxPyModuleId, "FilterLowScoreCalls" );
    exportModule<FilterDiagonalLineCalls>( rxPyModuleId, "FilterDiagonalLineCalls" );
    exportModule<FilterFuzzyCalls>( rxPyModuleId, "FilterFuzzyCalls" );
    exportModule<ComputeCallAmbiguity>( rxPyModuleId, "ComputeCallAmbiguity" );
} // function
#endif