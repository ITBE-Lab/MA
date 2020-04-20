#include "ma/module/compare_alignments.h"
#include "ms/util/pybind11.h"
#include <pybind11/stl.h>

using namespace libMA;
using namespace libMS;

#ifdef WITH_PYTHON

void exportCompareAlignments( libMS::SubmoduleOrganizer& xOrganizer )
{
    exportModule<AlignmentToSeeds>( xOrganizer, "AlignmentToSeeds" );
    exportModule<CompareSeedSets>( xOrganizer, "CompareSeedSets" );
    exportModule<CollectSeedSetComps>( xOrganizer, "CollectSeedSetComps", []( auto&& x ) {
        x.def_readwrite( "collection", &CollectSeedSetComps::pCollection );
    } );

    py::class_<SeedSetComp, libMS::Container, std::shared_ptr<SeedSetComp>>( xOrganizer.container( ), "SeedSetComp" )
        .def( py::init<>( ) )
        .def_readwrite( "nt_ground_truth", &SeedSetComp::uiNtGroundTruth )
        .def_readwrite( "nt_overlap", &SeedSetComp::uiNtOverlap )
        .def_readwrite( "nt_data", &SeedSetComp::uiNtData )
        .def_readwrite( "amount_ground_truth", &SeedSetComp::uiAmountGroundTruth )
        .def_readwrite( "amount_overlap", &SeedSetComp::uiAmountOverlap )
        .def_readwrite( "amount_90_percent_overlap", &SeedSetComp::uiAmount90PercentOverlap )
        .def_readwrite( "amount_data", &SeedSetComp::uiAmountData )
        .def_readwrite( "seeds_found", &SeedSetComp::xSeedsFound )
        .def( "add_ground_truth", &SeedSetComp::addGroundTruth )
        .def( "merge", &SeedSetComp::merge );
} // function

#endif