/**
 * @file seed.cpp
 * @author Markus Schmidt
 */
#include "container/seed.h"
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

#ifdef WITH_PYTHON

void exportSeed( py::module& rxPyModuleId )
{
    // export the Seed class
    py::class_<Seed>( rxPyModuleId, "Seed" )
        .def( py::init<nucSeqIndex, nucSeqIndex, nucSeqIndex, bool>() )
        .def_readwrite( "start", &Seed::iStart )
        .def_readwrite( "size", &Seed::iSize )
        .def_readwrite( "delta", &Seed::uiDelta )
        .def_readwrite( "start_ref", &Seed::uiPosOnReference )
        .def_readwrite( "on_forward_strand", &Seed::bOnForwStrand )
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
        .def( "extractStrand", &Seeds::extractStrand)
        .def( "extend", &Seeds::append)
        .def( "sort_by_ref_pos", &Seeds::sortByRefPos);

    // tell boost python that pointers of these classes can be converted implicitly
    py::implicitly_convertible<Seeds, Container>( );
} // function
#endif