#include "ma/container/alignment.h"
#include "ma/container/pack.h"
#include "ma/container/seed.h"
#include "ma/module/harmonization.h"
#include "ma/module/sam_reader.h"

#pragma once

namespace libMA
{

class AlignmentToSeeds : public libMS::Module<Seeds, false, Alignment>
{
  public:
    /**
     * @brief creates a new AlignmentToSeeds module.
     */
    AlignmentToSeeds( const ParameterSetManager& rParameters )
    {} // constructor

    std::shared_ptr<Seeds> DLL_PORT( MA ) execute( std::shared_ptr<Alignment> pAlignment )
    {
        return pAlignment->toSeeds( );
    } // method

}; // class AlignmentToSeeds


/**
 * @brief compares two seed sets
 * @details
 * Expects both seed sets to be sorted & lumped by SeedLumping.
 * Iterates over the sorted & lumped seed sets and computes the total overlap between seeds (as amount of nt).
 * A lumped seed set can never have two seeds overlapping within itself therefore no check is done for that.
 * This can deal with several seeds of one set overlapping a single seed of the other set though.
 */
class CompareSeedSets : public libMS::Module<SeedSetComp, false, Seeds, Seeds, SeedSetComp>
{
  public:
    /**
     * @brief returns wether xA is considered smaller than xB
     */
    static bool isSmaller( const Seed& xA, const Seed& xB )
    {
        if( xA.bOnForwStrand != xB.bOnForwStrand )
            return xA.bOnForwStrand;
        if( SeedLumping::getDelta( xA ) != SeedLumping::getDelta( xB ) )
            return SeedLumping::getDelta( xA ) < SeedLumping::getDelta( xB );
        return xA.end( ) <= xB.start( );
    }

    /**
     * @brief creates a new CompareSeedSets module.
     */
    CompareSeedSets( const ParameterSetManager& rParameters )
    {} // constructor

    std::shared_ptr<SeedSetComp> DLL_PORT( MA )
        execute( std::shared_ptr<Seeds> pGroundTruth, std::shared_ptr<Seeds> pData, std::shared_ptr<SeedSetComp> pComp )
    {
        auto xGroundTruthIt = pGroundTruth->begin( );
        auto xDataIt = pData->begin( );

        pComp->uiAmountData += pData->size( );

        while( xGroundTruthIt != pGroundTruth->end( ) && xDataIt != pData->end( ) )
        {
            if( isSmaller( *xDataIt, *xGroundTruthIt ) )
                pComp->addData( xDataIt );
            else if( isSmaller( *xGroundTruthIt, *xDataIt ) )
                ++xGroundTruthIt;
            else
            {
                pComp->addOverlap( xDataIt, xGroundTruthIt );
                if( xDataIt->end( ) < xGroundTruthIt->end( ) )
                    pComp->addData( xDataIt );
                else
                    ++xGroundTruthIt;
            } // else
        } // while
        // remaining seeds in data
        while( xDataIt != pData->end( ) )
            pComp->addData( xDataIt );

        return pComp;
    } // method
}; // class CompareSeedSets


class CollectSeedSetComps : public libMS::Module<libMS::Container, false, SeedSetComp>
{
  public:
    std::shared_ptr<SeedSetComp> pCollection;
    std::mutex xMutex;

    /**
     * @brief creates a new CollectSeedSetComps module.
     */
    CollectSeedSetComps( const ParameterSetManager& rParameters ) : pCollection( std::make_shared<SeedSetComp>( ) )
    {} // constructor

    std::shared_ptr<libMS::Container> DLL_PORT( MA ) execute( std::shared_ptr<SeedSetComp> pNew )
    {
        std::lock_guard<std::mutex> xGuard( xMutex );
        pCollection->merge( *pNew );
        return std::make_shared<libMS::Container>( );
    } // method

}; // class CollectSeedSetComps

} // namespace libMA

#ifdef WITH_PYTHON
void exportCompareAlignments( libMS::SubmoduleOrganizer& xOrganizer );
#endif