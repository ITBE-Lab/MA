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
        rxPyModuleId, "SvCallSink", []( auto&& x ) {
            x.def_readwrite( "run_id", &SvCallSink::iRunId );
        } );


    exportModule<CompleteBipartiteSubgraphSweep, std::shared_ptr<SV_DB>, std::shared_ptr<Pack>, int64_t>(
        rxPyModuleId, "CompleteBipartiteSubgraphSweep" );

    exportModule<ExactCompleteBipartiteSubgraphSweep, std::shared_ptr<SV_DB>, std::shared_ptr<Pack>, int64_t>(
        rxPyModuleId, "ExactCompleteBipartiteSubgraphSweep" );
} // function
#endif