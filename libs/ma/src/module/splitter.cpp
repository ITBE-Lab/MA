/**
 * @file splitter.cpp
 * @author Markus Schmidt
 */
#include "ma/module/splitter.h"
#include "ma/container/alignment.h"
#include "ma/container/nucSeq.h"
#include "ma/container/soc.h"
#include "ma/module/filter_seeds_by_area.h"
#include "ms/module/splitter.h"

using namespace libMS;
using namespace libMA;


#ifdef WITH_PYTHON
#include "ms/util/pybind11.h"
#include "pybind11/pybind11.h"

void exportSplitter( libMS::SubmoduleOrganizer& xOrganizer )
{
    // export the Lock class
    exportModule<Lock<libMS::Container>>( xOrganizer, "Lock" );

    // export the UnLock class
    exportModule<UnLock<libMS::Container>, std::shared_ptr<BasePledge>>( xOrganizer, "UnLock" );

    // export the Splitter<NucSeq> class
    exportModule<Splitter<NucSeq>>( xOrganizer, "NucSeqSplitter" );
    // export the StaticSplitter<NucSeq> class
    exportModule<StaticSplitter<NucSeq>, std::shared_ptr<ContainerVector<std::shared_ptr<NucSeq>>>>(
        xOrganizer, "StaticNucSeqSplitter" );

    // exported in fileReader.cpp
    // py::bind_vector_ext<ContainerVector<std::shared_ptr<NucSeq>>, libMS::Container,
    //                    std::shared_ptr<ContainerVector<std::shared_ptr<NucSeq>>>>( xOrganizer,
    //                                                                                "NucSeqContainerVector" );

    // export the TupleGet class
    exportModule<TupleGet<ContainerVector<std::shared_ptr<NucSeq>>, 0>>( xOrganizer, "GetFirstQuery" );
    exportModule<TupleGet<ContainerVector<std::shared_ptr<NucSeq>>, 1>>( xOrganizer, "GetSecondQuery" );

    py::bind_vector<std::vector<std::tuple<std::shared_ptr<NucSeq>, std::shared_ptr<SoCPriorityQueue>>>>(
        xOrganizer._util(), "SoCPriorityQueueVector", "docstr" );

    exportModule<Collector<NucSeq, SoCPriorityQueue>>( xOrganizer, "NucSeqSoCCollector", []( auto&& x ) {
        x.def_readwrite( "collection", &Collector<NucSeq, SoCPriorityQueue>::vCollection );
    } );

    py::bind_vector<std::vector<std::tuple<std::shared_ptr<NucSeq>, std::shared_ptr<SegmentVector>>>>(
        xOrganizer._util(), "NucSeqSegmentVector", "docstr" );
    exportModule<Collector<NucSeq, SegmentVector>>( xOrganizer, "NucSeqSegmentCollector", []( auto&& x ) {
        x.def_readwrite( "collection", &Collector<NucSeq, SegmentVector>::vCollection );
    } );

    exportModule<Collector<NucSeq, NucSeq, SoCPriorityQueue, SoCPriorityQueue>>(
        xOrganizer, "PairedNucSeqSoCCollector", []( auto&& x ) {
            x.def_readwrite( "collection",
                             &Collector<NucSeq, NucSeq, SoCPriorityQueue, SoCPriorityQueue>::vCollection );
        } );

    exportModule<Collector<NucSeq, ContainerVector<std::shared_ptr<Seeds>>>>(
        xOrganizer, "NucSeqSeedsCollector", []( auto&& x ) {
            x.def_readwrite( "collection", &Collector<NucSeq, ContainerVector<std::shared_ptr<Seeds>>>::vCollection );
        } );

    exportModule<Join<libMS::Container, libMS::Container>>( xOrganizer, "ContainerJoin" );

    exportModule<Collector<NucSeq, ContainerVector<std::shared_ptr<Alignment>>, Pack>>(
        xOrganizer, "AlignmentCollector", []( auto&& x ) {
            x.def_readwrite( "collection",
                             &Collector<NucSeq, ContainerVector<std::shared_ptr<Alignment>>, Pack>::vCollection );
        } );

    exportModule<FilterSeedsByArea, nucSeqIndex, nucSeqIndex>( xOrganizer, "FilterSeedsByArea" );
    exportModule<VectorCollector<Seeds>>( xOrganizer, "SeedsCollector", []( auto&& x ) {
        x.def_readwrite( "collection", &VectorCollector<Seeds>::pCollection );
    } );
} // function
#endif
