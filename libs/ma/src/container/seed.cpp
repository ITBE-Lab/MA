/**
 * @file seed.cpp
 * @author Markus Schmidt
 */
#include "container/seed.h"
#include "container/pack.h"
#include "util/default_parameters.h"
#include "util/pybind11.h"
using namespace libMA;

#ifdef _MSC_VER
// getting ambiguous abs otherwise
#include <cstdlib>
#endif

using namespace libMA::defaults;
extern int libMA::defaults::iGap;
extern int libMA::defaults::iExtend;
extern int libMA::defaults::iMatch;
extern int libMA::defaults::iMissMatch;

void Seeds::confirmSeedPositions( std::shared_ptr<NucSeq> pQuery, std::shared_ptr<Pack> pRef, bool bIsMaxExtended )
{
    static const char chars[ 5 ] = {'A', 'C', 'G', 'T', 'N'};
    size_t uiI = 0;
    for( auto& rSeed : vContent )
    {
        std::shared_ptr<NucSeq> pRefSec = std::make_shared<NucSeq>( );
        if( pRef->bridgingPositions( rSeed.start_ref( ), rSeed.end_ref( ) ) )
            continue;
        pRef->vExtractSubsectionN( rSeed.start_ref( ), rSeed.end_ref( ), *pRefSec );
        // pRef->vExtractSubsectionN( pRef->uiPositionToReverseStrand( rSeed.end_ref( ) ),
        //                           pRef->uiPositionToReverseStrand( rSeed.start_ref( ) ),
        //                           *pRefSec );
        uint8_t uiBefore = 4;
        if( rSeed.start_ref( ) > 0 )
            uiBefore = pRef->isHole( pRef->iAbsolutePosition( rSeed.start_ref( ) - 1 ) )
                           ? 4
                           : pRef->vExtract( rSeed.start_ref( ) - 1 );

        uint8_t uiAfter = 4;
        if( rSeed.end_ref( ) < pRef->uiUnpackedSizeForwardStrand * 2 )
            uiAfter =
                pRef->isHole( pRef->iAbsolutePosition( rSeed.end_ref( ) ) ) ? 4 : pRef->vExtract( rSeed.end_ref( ) );

        bool bFailed = false;
        for( size_t uiI = 0; uiI < rSeed.size( ); uiI++ )
            if( pQuery->pxSequenceRef[ rSeed.start( ) + uiI ] != pRefSec->pxSequenceRef[ uiI ] )
                bFailed = true;
        if( bIsMaxExtended && rSeed.start( ) > 0 && rSeed.start_ref( ) > 0 &&
            pQuery->pxSequenceRef[ rSeed.start( ) - 1 ] == uiBefore )
            bFailed = true;
        if( bIsMaxExtended && rSeed.end( ) < pQuery->length( ) &&
            rSeed.end_ref( ) < pRef->uiUnpackedSizeForwardStrand * 2 &&
            pQuery->pxSequenceRef[ rSeed.end( ) ] == uiAfter )
            bFailed = true;
        if( bFailed )
        {
            std::cout << "bIsMaxExtended = " << ( bIsMaxExtended ? "true" : "false" ) << std::endl;
            std::cout << "Seed: (number " << uiI << ")" << std::endl;
            std::cout << "(q,r,l): " << rSeed.start( ) << ", " << rSeed.start_ref( ) << ", " << rSeed.size( )
                      << std::endl;
            std::cout << "q_len " << pQuery->length( ) << " ref_len " << pRef->uiUnpackedSizeForwardStrand << std::endl;
            std::cout << "Query:" << std::endl;
            if( rSeed.start( ) > 1 )
                std::cout << chars[ pQuery->pxSequenceRef[ rSeed.start( ) - 2 ] ];
            else
                std::cout << " ";
            if( rSeed.start( ) > 0 )
                std::cout << chars[ pQuery->pxSequenceRef[ rSeed.start( ) - 1 ] ] << " ";
            else
                std::cout << "  ";
            for( nucSeqIndex uiI = rSeed.start( ); uiI < rSeed.end( ); uiI++ )
                std::cout << chars[ pQuery->pxSequenceRef[ uiI ] ];
            if( rSeed.end( ) < pQuery->length( ) )
                std::cout << " " << chars[ pQuery->pxSequenceRef[ rSeed.end( ) ] ];
            if( rSeed.end( ) + 1 < pQuery->length( ) )
                std::cout << chars[ pQuery->pxSequenceRef[ rSeed.end( ) + 1 ] ];
            std::cout << std::endl;
            std::cout << "Ref:" << std::endl;
            if( rSeed.start_ref( ) > 1 )
                std::cout << chars[ pRef->isHole( pRef->iAbsolutePosition( rSeed.start_ref( ) - 2 ) )
                                        ? 4
                                        : pRef->vExtract( rSeed.start_ref( ) - 2 ) ];
            else
                std::cout << " ";
            if( rSeed.start_ref( ) > 0 )
                std::cout << chars[ uiBefore ] << " ";
            else
                std::cout << "  ";
            for( nucSeqIndex uiI = 0; uiI < rSeed.size( ); uiI++ )
                std::cout << chars[ pRefSec->pxSequenceRef[ uiI ] ];
            if( rSeed.end_ref( ) < pRef->uiUnpackedSizeForwardStrand * 2 )
                std::cout << " " << chars[ uiAfter ];
            if( rSeed.end_ref( ) + 1 < pRef->uiUnpackedSizeForwardStrand * 2 )
                std::cout << chars[ pRef->isHole( pRef->iAbsolutePosition( rSeed.end_ref( ) + 1 ) )
                                        ? 4
                                        : pRef->vExtract( rSeed.end_ref( ) + 1 ) ];
            std::cout << std::endl;
            if( rSeed.start_ref( ) >= 2 && rSeed.end_ref( ) + 2 <= pRef->uiUnpackedSizeForwardStrand * 2 )
            {
                std::cout << "ref +-2" << std::endl;
                std::shared_ptr<NucSeq> pRefSec2 = std::make_shared<NucSeq>( );
                pRef->vExtractSubsectionN( rSeed.start_ref( ) - 2, rSeed.end_ref( ) + 2, *pRefSec2 );
                std::cout << pRefSec2->toString( ) << std::endl;
            } // if

            throw std::runtime_error( "computed wrong seed!" );
        } // if
        uiI++;
    } // for
} // method

#ifdef WITH_PYTHON

#ifdef BOOST_PYTHON
void exportSeed( )
{
    exportInterval<nucSeqIndex>( );
    // export the Seed class
    boost::python::class_<Seed, boost::python::bases<Interval<nucSeqIndex>>>( "Seed" )
        .def_readwrite( "start_ref", &Seed::uiPosOnReference )
        .def( "__eq__", &Seed::operator==);


    // export the Seed class
    boost::python::class_<AlignmentStatistics>( "AlignmentStatistics", boost::python::init<>( ) )
        .def_readwrite( "index_of_strip", &AlignmentStatistics::index_of_strip )
        .def_readwrite( "num_seeds_in_strip", &AlignmentStatistics::num_seeds_in_strip )
        .def_readwrite( "anchor_size", &AlignmentStatistics::anchor_size )
        .def_readwrite( "anchor_ambiguity", &AlignmentStatistics::anchor_ambiguity )
        .def_readwrite( "name", &AlignmentStatistics::sName )
        .def_readwrite( "initial_q_beg", &AlignmentStatistics::uiInitialQueryBegin )
        .def_readwrite( "initial_r_beg", &AlignmentStatistics::uiInitialRefBegin )
        .def_readwrite( "initial_q_end", &AlignmentStatistics::uiInitialQueryEnd )
        .def_readwrite( "initial_r_end", &AlignmentStatistics::uiInitialRefEnd );

    // export the Seeds class
    boost::python::
        class_<Seeds, boost::python::bases<Container>, boost::python::bases<std::list<Seed>>, std::shared_ptr<Seeds>>(
            "Seeds" )
            .def( boost::python::init<std::shared_ptr<Seeds>>( ) )
            .def( boost::python::vector_indexing_suite<Seeds,
                                                       /*
                                                        *    true = noproxy this means that the content of
                                                        * the vector is already exposed by boost python.
                                                        *    if this is kept as false, Container would be
                                                        * exposed a second time. the two Containers would
                                                        * be different and not inter castable.
                                                        */
                                                       true>( ) );

    // make vectors of container-pointers a thing
    IterableConverter( ).from_python<Seeds>( );

    // tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible<std::shared_ptr<Seeds>, std::shared_ptr<Container>>( );
} // function
#else
void exportSeed( py::module& rxPyModuleId )
{
    // export the Seed class
    py::class_<Seed>( rxPyModuleId, "Seed" )
        .def_readwrite( "start_ref", &Seed::uiPosOnReference )
        .def_readwrite( "start", &Seed::iStart )
        .def_readwrite( "size", &Seed::iSize )
        .def( "__eq__", &Seed::operator==);

    py::class_<AlignmentStatistics>( rxPyModuleId, "AlignmentStatistics" )
        .def( py::init<>( ) )
        .def_readwrite( "index_of_strip", &AlignmentStatistics::index_of_strip )
        .def_readwrite( "num_seeds_in_strip", &AlignmentStatistics::num_seeds_in_strip )
        .def_readwrite( "anchor_size", &AlignmentStatistics::anchor_size )
        .def_readwrite( "anchor_ambiguity", &AlignmentStatistics::anchor_ambiguity )
        .def_readwrite( "name", &AlignmentStatistics::sName )
        .def_readwrite( "initial_q_beg", &AlignmentStatistics::uiInitialQueryBegin )
        .def_readwrite( "initial_r_beg", &AlignmentStatistics::uiInitialRefBegin )
        .def_readwrite( "initial_q_end", &AlignmentStatistics::uiInitialQueryEnd )
        .def_readwrite( "initial_r_end", &AlignmentStatistics::uiInitialRefEnd );

    // export the Seeds class
    py::bind_vector_ext<Seeds, Container, std::shared_ptr<Seeds>>( rxPyModuleId, "Seeds", "docstr" )
        .def( py::init<std::shared_ptr<Seeds>>( ) )
        .def( py::init<>( ) )
        .def( "compare_seed_sets", &Seeds::compareSeedSets )
        .def( "split_seed_sets", &Seeds::splitSeedSets )
        .def( "confirm_seed_positions", &Seeds::confirmSeedPositions )
        .def( "sort_by_ref_pos", &Seeds::sortByRefPos )
        .def( "sort_by_q_pos", &Seeds::sortByQPos );

    // tell boost python that pointers of these classes can be converted implicitly
    py::implicitly_convertible<Seeds, Container>( );
} // function
#endif
#endif