/**
 * @file sweepSvJumps.cpp
 * @author Markus Schmidt
 */
#include "module/sweepSvJumps.h"
#include "util/pybind11.h"

using namespace libMA;


#ifdef WITH_PYTHON
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
    exportModule<SvCallSink, std::shared_ptr<SV_DB>, std::string, std::string, int64_t>(
        rxPyModuleId, "SvCallSink", []( auto&& x ) { x.def_readwrite( "run_id", &SvCallSink::iRunId ); } );
    exportModule<BufferedSvCallSink, std::shared_ptr<SV_DB>, int64_t>(
        rxPyModuleId, "BufferedSvCallSink", []( auto&& x ) { x.def( "commit", &BufferedSvCallSink::commit ); } );


    exportModule<CompleteBipartiteSubgraphSweep, std::shared_ptr<SV_DB>, std::shared_ptr<Pack>, int64_t, int64_t>(
        rxPyModuleId, "CompleteBipartiteSubgraphSweep", []( auto&& x ) {
            x.def_readonly( "time_init", &CompleteBipartiteSubgraphSweep::dInit )
                .def_readonly( "time_complete_while", &CompleteBipartiteSubgraphSweep::dOuterWhile )
                .def_readonly( "time_inner_while", &CompleteBipartiteSubgraphSweep::dInnerWhile );
        } );

    exportModule<ExactCompleteBipartiteSubgraphSweep, std::shared_ptr<SV_DB>, std::shared_ptr<Pack>, int64_t>(
        rxPyModuleId, "ExactCompleteBipartiteSubgraphSweep" );
    exportModule<FilterLowSupportShortCalls>( rxPyModuleId, "FilterLowSupportShortCalls" );
    exportModule<FilterDiagonalLineCalls>( rxPyModuleId, "FilterDiagonalLineCalls" );
    exportModule<FilterLowCoverageCalls, std::shared_ptr<SV_DB>, std::shared_ptr<Pack>, int64_t>(
        rxPyModuleId, "FilterLowCoverageCalls" );
    exportModule<FilterFuzzyCalls>( rxPyModuleId, "FilterFuzzyCalls" );
    exportModule<ComputeCallAmbiguity>( rxPyModuleId, "ComputeCallAmbiguity" );
} // function
#endif