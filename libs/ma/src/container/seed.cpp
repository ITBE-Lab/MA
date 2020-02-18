/**
 * @file seed.cpp
 * @author Markus Schmidt
 */
#include "container/seed.h"
#include "container/pack.h"
#include "util/pybind11.h"
using namespace libMA;

#ifdef _MSC_VER
// getting ambiguous abs otherwise
#include <cstdlib>
#endif

void Seeds::confirmSeedPositions( std::shared_ptr<NucSeq> pQuery, std::shared_ptr<Pack> pRef, bool bIsMaxExtended )
{
    return; //@todo fix 'n' and other symbols in reference...
    static const char chars[ 5 ] = {'A', 'C', 'G', 'T', 'N'};
    for( auto& rSeed : vContent )
    {
        std::shared_ptr<NucSeq> pRefSec = std::make_shared<NucSeq>( );
        if( rSeed.bOnForwStrand )
        {
            if( pRef->bridgingPositions( rSeed.start_ref( ), rSeed.end_ref( ) ) )
                continue;
            pRef->vExtractSubsectionN( rSeed.start_ref( ), rSeed.end_ref( ), *pRefSec );
        } // if
        else
        {
            if( pRef->bridgingPositions( rSeed.start_ref( ) - rSeed.size( ) + 1, rSeed.start_ref( ) + 1 ) )
                continue;
            pRef->vExtractSubsectionN( rSeed.start_ref( ) - rSeed.size( ) + 1, rSeed.start_ref( ) + 1, *pRefSec );
            pRefSec->vReverseAll( );
            pRefSec->vSwitchAllBasePairsToComplement( );
        }
        // pRef->vExtractSubsectionN( pRef->uiPositionToReverseStrand( rSeed.end_ref( ) ),
        //                           pRef->uiPositionToReverseStrand( rSeed.start_ref( ) ),
        //                           *pRefSec );
        uint8_t uiBefore = 4;
        if( rSeed.start_ref( ) > 0 && rSeed.bOnForwStrand )
            uiBefore = pRef->isHole( rSeed.start_ref( ) - 1 ) ? 4 : pRef->getNucleotideOnPos( rSeed.start_ref( ) - 1 );

        auto uiX = rSeed.start_ref( ) + 1;
        if( uiX < pRef->uiUnpackedSizeForwardStrand && !rSeed.bOnForwStrand )
            uiBefore = pRef->isHole( uiX ) ? 4 : 3 - pRef->getNucleotideOnPos( uiX );

        uint8_t uiAfter = 4;
        if( rSeed.end_ref( ) < pRef->uiUnpackedSizeForwardStrand && rSeed.bOnForwStrand )
            uiAfter = pRef->isHole( rSeed.end_ref( ) ) ? 4 : pRef->getNucleotideOnPos( rSeed.end_ref( ) );
        uiX = rSeed.start_ref( ) - rSeed.size( ) + 1;
        if( uiX > 0 && !rSeed.bOnForwStrand )
            uiAfter = pRef->isHole( uiX - 1 ) ? 4 : 3 - pRef->getNucleotideOnPos( uiX - 1 );

        bool bFailed = false;
        for( size_t uiI = 0; uiI < rSeed.size( ); uiI++ )
            if( pQuery->pxSequenceRef[ rSeed.start( ) + uiI ] != pRefSec->pxSequenceRef[ uiI ] )
                bFailed = true;
        if( bIsMaxExtended && rSeed.start( ) > 0 && rSeed.start_ref( ) > 0 &&
            pQuery->pxSequenceRef[ rSeed.start( ) - 1 ] == uiBefore )
            bFailed = true;
        if( bIsMaxExtended && rSeed.end( ) < pQuery->length( ) &&
            rSeed.end_ref( ) < pRef->uiUnpackedSizeForwardStrand && pQuery->pxSequenceRef[ rSeed.end( ) ] == uiAfter )
            bFailed = true;
        if( bFailed )
        {
            std::cout << "Query:" << std::endl;
            if( rSeed.start( ) > 0 )
                std::cout << chars[ pQuery->pxSequenceRef[ rSeed.start( ) - 1 ] ] << " ";
            else
                std::cout << "  ";
            for( nucSeqIndex uiI = rSeed.start( ); uiI < rSeed.end( ); uiI++ )
                std::cout << chars[ pQuery->pxSequenceRef[ uiI ] ];
            if( rSeed.end( ) < pQuery->length( ) )
                std::cout << " " << chars[ pQuery->pxSequenceRef[ rSeed.end( ) ] ];
            std::cout << std::endl;
            std::cout << "Ref:" << std::endl;
            if( rSeed.start_ref( ) > 0 )
                std::cout << chars[ uiBefore ] << " ";
            else
                std::cout << "  ";
            for( nucSeqIndex uiI = 0; uiI < rSeed.size( ); uiI++ )
                std::cout << chars[ pRefSec->pxSequenceRef[ uiI ] ];
            if( rSeed.end_ref( ) < pRef->uiUnpackedSizeForwardStrand )
                std::cout << " " << chars[ uiAfter ];
            std::cout << std::endl;

            throw std::runtime_error( "computed wrong seed!" );
        } // if
    } // for
} // method

#ifdef WITH_PYTHON

void exportSeed( py::module& rxPyModuleId )
{
    // export the Seed class
    py::class_<Seed>( rxPyModuleId, "Seed" )
        .def( py::init<nucSeqIndex, nucSeqIndex, nucSeqIndex, bool>( ) )
        .def_readwrite( "start", &Seed::iStart )
        .def_readwrite( "size", &Seed::iSize )
        .def_readwrite( "delta", &Seed::uiDelta )
        .def_readwrite( "start_ref", &Seed::uiPosOnReference )
        .def_readwrite( "on_forward_strand", &Seed::bOnForwStrand )
#if DEBUG_LEVEL > 0
        .def_readwrite( "id", &Seed::uiId )
#endif
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
        .def( "extractStrand", &Seeds::extractStrand )
        .def( "extend", &Seeds::append )
        .def( "compare_seed_sets", &Seeds::compareSeedSets )
        .def( "split_seed_sets", &Seeds::splitSeedSets )
        .def( "confirm_seed_positions", &Seeds::confirmSeedPositions )
        .def( "sort_by_ref_pos", &Seeds::sortByRefPos )
        .def( "sort_by_q_pos", &Seeds::sortByQPos );

    // tell boost python that pointers of these classes can be converted implicitly
    py::implicitly_convertible<Seeds, Container>( );
} // function
#endif