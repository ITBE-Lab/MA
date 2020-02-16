/**
 * @file splitter.cpp
 * @author Markus Schmidt
 */
#include "module/splitter.h"

using namespace libMA;


#ifdef WITH_PYTHON
#include "util/pybind11.h"

void exportSplitter( py::module& rxPyModuleId )
{
    // export the Lock class
    exportModule<Lock<Container>>( rxPyModuleId, "Lock" );

    // export the UnLock class
    exportModule<UnLock<Container>, std::shared_ptr<BasePledge>>( rxPyModuleId, "UnLock" );

    // export the Splitter<NucSeq> class
    exportModule<Splitter<NucSeq>>( rxPyModuleId, "NucSeqSplitter" );
    // export the StaticSplitter<NucSeq> class
    exportModule<StaticSplitter<NucSeq>, std::shared_ptr<ContainerVector<std::shared_ptr<NucSeq>>>>(
        rxPyModuleId, "StaticNucSeqSplitter" );

    // exported in fileReader.cpp
    // py::bind_vector_ext<ContainerVector<std::shared_ptr<NucSeq>>, Container,
    //                    std::shared_ptr<ContainerVector<std::shared_ptr<NucSeq>>>>( rxPyModuleId,
    //                                                                                "NucSeqContainerVector" );

    // export the TupleGet class
    exportModule<TupleGet<ContainerVector<std::shared_ptr<NucSeq>>, 0>>( rxPyModuleId, "GetFirstQuery" );
    exportModule<TupleGet<ContainerVector<std::shared_ptr<NucSeq>>, 1>>( rxPyModuleId, "GetSecondQuery" );

    py::bind_vector<std::vector<std::tuple<std::shared_ptr<NucSeq>, std::shared_ptr<SoCPriorityQueue>>>>(
        rxPyModuleId, "SoCPriorityQueueVector", "docstr" );

    exportModule<Collector<NucSeq, SoCPriorityQueue>>( rxPyModuleId, "NucSeqSoCCollector", []( auto&& x ) {
        x.def_readwrite( "collection", &Collector<NucSeq, SoCPriorityQueue>::vCollection );
    } );

    py::bind_vector<std::vector<std::tuple<std::shared_ptr<NucSeq>, std::shared_ptr<SegmentVector>>>>(
        rxPyModuleId, "NucSeqSegmentVector", "docstr" );
    exportModule<Collector<NucSeq, SegmentVector>>( rxPyModuleId, "NucSeqSegmentCollector", []( auto&& x ) {
        x.def_readwrite( "collection", &Collector<NucSeq, SegmentVector>::vCollection );
    } );

    exportModule<Collector<NucSeq, NucSeq, SoCPriorityQueue, SoCPriorityQueue>>(
        rxPyModuleId, "PairedNucSeqSoCCollector", []( auto&& x ) {
            x.def_readwrite( "collection",
                             &Collector<NucSeq, NucSeq, SoCPriorityQueue, SoCPriorityQueue>::vCollection );
        } );

    exportModule<Collector<NucSeq, ContainerVector<std::shared_ptr<Seeds>>>>(
        rxPyModuleId, "NucSeqSeedsCollector", []( auto&& x ) {
            x.def_readwrite( "collection", &Collector<NucSeq, ContainerVector<std::shared_ptr<Seeds>>>::vCollection );
        } );

    exportModule<Join<Container, Container>>( rxPyModuleId, "ContainerJoin" );

    exportModule<Collector<NucSeq, ContainerVector<std::shared_ptr<Alignment>>, Pack>>(
        rxPyModuleId, "AlignmentCollector", []( auto&& x ) {
            x.def_readwrite( "collection",
                             &Collector<NucSeq, ContainerVector<std::shared_ptr<Alignment>>, Pack>::vCollection );
        } );

    exportModule<FilterSeedsByArea, nucSeqIndex, nucSeqIndex>( rxPyModuleId, "FilterSeedsByArea" );
    exportModule<VectorCollector<Seeds>>( rxPyModuleId, "SeedsCollector", []( auto&& x ) {
        x.def_readwrite( "collection", &VectorCollector<Seeds>::pCollection );
    } );
} // function
#endif
