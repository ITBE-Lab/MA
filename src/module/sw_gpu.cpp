/**
 * @file execOnVector.cpp
 * @author Markus Schmidt
 */
#include "module/sw_gpu.h"
#include "container/nucSeq.h"
#include "module/module.h"
using namespace libMA;

#ifdef WITH_GPU_SW
std::vector<GPUReturn> testGPUSW( std::shared_ptr<ContainerVector> pQueries,
                                  std::shared_ptr<NucSeq> b, unsigned int uiGpuId )
{
    std::vector<std::vector<char>> vQueries;
    for ( const auto& pQuery : *pQueries )
    {
        const auto& a = std::dynamic_pointer_cast<NucSeq>( pQuery ); // dc
        std::vector<char> vQuery( a->pGetSequenceRef( ), a->pGetSequenceRef( ) + a->length( ) );
        vQueries.push_back( vQuery );
    } // for
    std::vector<char> vRef( b->pGetSequenceRef( ), b->pGetSequenceRef( ) + b->length( ) );

    auto vResults = cudaAlign( vRef, vQueries, uiGpuId );

    return vResults;
} // function
#endif

void exportSW_GPU( )
{
#ifdef WITH_GPU_SW
    boost::python::def( "testGPUSW", &testGPUSW );

    // export the Tail class
    boost::python::class_<GPUReturn>( "GPUReturn" )
        .def_readwrite( "vMaxPos", &GPUReturn::vMaxPos )
        .def_readwrite( "iMaxScore", &GPUReturn::iMaxScore );

    boost::python::class_<std::vector<GPUReturn>>( "GPUReturnList" )
        .def(
            boost::python::vector_indexing_suite<std::vector<GPUReturn>,
                                                 /*
                                                  *    true = noproxy this means that the content of
                                                  * the vector is already exposed by boost python.
                                                  *    if this is kept as false, Container would be
                                                  * exposed a second time. the two Containers would
                                                  * be different and not inter castable.
                                                  */
                                                 true>( ) );

    boost::python::class_<std::vector<size_t>>( "GPUReturnMaxPosList" )
        .def(
            boost::python::vector_indexing_suite<std::vector<size_t>,
                                                 /*
                                                  *    true = noproxy this means that the content of
                                                  * the vector is already exposed by boost python.
                                                  *    if this is kept as false, Container would be
                                                  * exposed a second time. the two Containers would
                                                  * be different and not inter castable.
                                                  */
                                                 true>( ) );
#endif
} // function