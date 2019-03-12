/**
 * @file hash_map_seeding.cpp
 * @author Markus Schmidt
 */
#include "module/hashMapSeeding.h"
#include <limits>
#include <unordered_map>

using namespace libMA;

using namespace libMA::defaults;


std::shared_ptr<Seeds> HashMapSeeding::execute( std::shared_ptr<NucSeq> pQ1, std::shared_ptr<NucSeq> pQ2 )
{
    std::unordered_multimap<std::string, size_t> xIndex;

    // insert seeds into index
    for( size_t uiI = 0; uiI < pQ2->length( ); uiI += 1 )
        xIndex.emplace( pQ2->fromTo( uiI, uiI + this->uiSeedSize ), uiI );

    auto pSeeds = std::make_shared<Seeds>( );
    for( size_t uiI = 0; uiI < pQ1->length( ); uiI += 1 )
    {
        auto tuiRange = xIndex.equal_range( pQ1->fromTo( uiI, uiI + this->uiSeedSize ) );
        for( auto xIt = tuiRange.first; xIt != tuiRange.second; ++xIt )
            pSeeds->emplace_back( uiI, this->uiSeedSize, xIt->second, true );
    } // for

    return pSeeds;
} // function

#ifdef WITH_PYTHON
void exportHashMapSeeding( py::module& rxPyModuleId )
{
    // export the HashMapSeeding class
    exportModule<HashMapSeeding>( rxPyModuleId, "HashMapSeeding" );
} // function
#endif
