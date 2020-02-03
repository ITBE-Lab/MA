/**
 * @file sweepSvJumps.cpp
 * @author Markus Schmidt
 */
#include "module/sweepSvJumps.h"
#include "util/pybind11.h"

using namespace libMA;


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
    exportModule<SvCallSink<DBCon>, std::shared_ptr<SV_Schema<DBCon>>, std::string, std::string, int64_t>(
        rxPyModuleId, "SvCallSink", []( auto&& x ) { x.def_readwrite( "run_id", &SvCallSink<DBCon>::iRunId ); } );
    exportModule<BufferedSvCallSink<DBCon>, std::shared_ptr<SvCallInserter<DBCon>>>(
        rxPyModuleId, "BufferedSvCallSink", []( auto&& x ) { x.def( "commit", &BufferedSvCallSink<DBCon>::commit ); } );


    exportModule<CompleteBipartiteSubgraphSweep<DBCon>,
                 std::shared_ptr<SV_Schema<DBCon>>,
                 std::shared_ptr<Pack>,
                 int64_t,
                 int64_t>( rxPyModuleId, "CompleteBipartiteSubgraphSweep", []( auto&& x ) {
        x.def_readonly( "time_init", &CompleteBipartiteSubgraphSweep<DBCon>::dInit )
            .def_readonly( "time_complete_while", &CompleteBipartiteSubgraphSweep<DBCon>::dOuterWhile )
            .def_readonly( "time_inner_while", &CompleteBipartiteSubgraphSweep<DBCon>::dInnerWhile );
    } );

    exportModule<ExactCompleteBipartiteSubgraphSweep<DBCon>,
                 std::shared_ptr<SV_Schema<DBCon>>,
                 std::shared_ptr<Pack>,
                 int64_t>( rxPyModuleId, "ExactCompleteBipartiteSubgraphSweep" );
    exportModule<FilterLowSupportShortCalls>( rxPyModuleId, "FilterLowSupportShortCalls" );
    exportModule<FilterLowScoreCalls>( rxPyModuleId, "FilterLowScoreCalls" );
    exportModule<FilterDiagonalLineCalls>( rxPyModuleId, "FilterDiagonalLineCalls" );
    exportModule<FilterFuzzyCalls>( rxPyModuleId, "FilterFuzzyCalls" );
    exportModule<ComputeCallAmbiguity>( rxPyModuleId, "ComputeCallAmbiguity" );
} // function
#endif