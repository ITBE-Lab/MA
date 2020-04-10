#include "ma/module/compare_alignments.h"
#include "ms/util/pybind11.h"

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
    py::class_<SeedSetComp, libMS::Container, std::shared_ptr<SeedSetComp>>( xOrganizer.container( ), "SeedSetComp" );
} // function

#endif