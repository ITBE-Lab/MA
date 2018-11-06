/**
 * @file splitter.cpp
 * @author Markus Schmidt
 */
#include "module/splitter.h"
#include "container/soc.h"

using namespace libMA;


#ifdef WITH_PYTHON
void exportSplitter( )
{
    // export the Lock class
    exportModule<Lock<Container>>( "Lock" );

    // export the UnLock class
    exportModule<UnLock<Container>, std::shared_ptr<BasePledge>>( "UnLock" );


    // export the TupleGet class
    exportModule<TupleGet<ContainerVector<std::shared_ptr<NucSeq>>, 0>>( "GetFirstQuery" );
    exportModule<TupleGet<ContainerVector<std::shared_ptr<NucSeq>>, 1>>( "GetSecondQuery" );

    typedef std::tuple<std::shared_ptr<NucSeq>, std::shared_ptr<SoCPriorityQueue>> TP_NUC_SOC;

    TupleToPython<TP_NUC_SOC>();
    boost::python::class_<std::vector<TP_NUC_SOC>, boost::noncopyable>( "SoCPriorityQueueVector",
                                                                        boost::python::no_init )
        /*
         * true = noproxy this means that the content of
         * the vector is already exposed by boost python.
         * if this is kept as false, Container would be
         * exposed a second time. the two Containers would
         * be different and not inter castable.
         */
        .def( boost::python::vector_indexing_suite<std::vector<TP_NUC_SOC>, true>( ) );

    exportModule<Collector<NucSeq, SoCPriorityQueue>>( "NucSeqSoCCollector", []( auto&& x ) {
        x.def_readwrite( "collection", &Collector<NucSeq, SoCPriorityQueue>::vCollection );
    } );
} // function
#endif
