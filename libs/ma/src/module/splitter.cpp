/**
 * @file splitter.cpp
 * @author Markus Schmidt
 */
#include "module/splitter.h"
#include "container/soc.h"

using namespace libMA;


#ifdef WITH_PYTHON
#include "util/pybind11.h"

#ifdef BOOST_PYTHON
void exportSplitter( )
{
    // export the Lock class
    exportModule<Lock<Container>>( "Lock" );

    // export the UnLock class
    exportModule<UnLock<Container>, std::shared_ptr<BasePledge>>( "UnLock" );


    // export the TupleGet class
    exportModule<TupleGet<ContainerVector<std::shared_ptr<NucSeq>>, 0>>( "GetFirstQuery" );
    exportModule<TupleGet<ContainerVector<std::shared_ptr<NucSeq>>, 1>>( "GetSecondQuery" );
} // function
#else
void exportSplitter( py::module& rxPyModuleId )
{
    // export the Lock class
    exportModule<Lock<Container>>( rxPyModuleId, "Lock" );

    // export the UnLock class
    exportModule<UnLock<Container>, std::shared_ptr<BasePledge>>( rxPyModuleId, "UnLock" );

    // export the Splitter<NucSeq> class
    exportModule<Splitter<NucSeq>>( rxPyModuleId, "NucSeqSplitter" );

    // exported in fileReader.cpp
    //py::bind_vector_ext<ContainerVector<std::shared_ptr<NucSeq>>, Container,
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
} // function
#endif
#endif
