#include "ma/container/alignment.h"
#include "ma/container/pack.h"
#include "ma/container/seed.h"
#include "ma/module/harmonization.h"

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


class SeedSetComp : public libMS::Container
{
  public:
    nucSeqIndex uiNtGroundTruth = 0;
    nucSeqIndex uiNtOverlap = 0;
    nucSeqIndex uiNtData = 0;

    size_t uiAmountGroundTruth = 0;
    size_t uiAmountOverlap = 0;
    size_t uiAmountData = 0;

    void addData( Seeds::iterator& xDataIt )
    {
        uiNtData += xDataIt->size( );
        xDataIt++;
    } // mehtod

    void addGroundTruth( Seeds::iterator& xGroundTruthIt )
    {
        uiNtGroundTruth += xGroundTruthIt->size( );
        xGroundTruthIt++;
    } // mehtod

    void addOverlap( Seeds::iterator& xDataIt, Seeds::iterator& xGroundTruthIt )
    {

        assert( std::min( xDataIt->end( ), xGroundTruthIt->end( ) ) >
                std::max( xDataIt->start( ), xGroundTruthIt->start( ) ) );
        uiNtOverlap += std::min( xDataIt->end( ), xGroundTruthIt->end( ) ) -
                       std::max( xDataIt->start( ), xGroundTruthIt->start( ) );
        uiAmountOverlap++;
    } // mehtod

    void merge( const SeedSetComp& xOther )
    {
        uiNtGroundTruth += xOther.uiNtGroundTruth;
        uiNtOverlap += xOther.uiNtOverlap;
        uiNtData += xOther.uiNtData;

        uiAmountGroundTruth += xOther.uiAmountGroundTruth;
        uiAmountOverlap += xOther.uiAmountOverlap;
        uiAmountData += xOther.uiAmountData;
    } // mehtod

}; // class SeedSetComp

/**
 * @brief compares two seed sets
 * @details
 * Expects both seed sets to be sorted & lumped by SeedLumping.
 * Iterates over the sorted & lumped seed sets and computes the total overlap between seeds (as amount of nt).
 * A lumped seed set can never have two seeds overlapping within itself therefore no check is done for that.
 * This can deal with several seeds of one set overlapping a single seed of the other set though.
 */
class CompareSeedSets : public libMS::Module<SeedSetComp, false, Seeds, Seeds>
{
  public:
    /**
     * @brief returns wether xA is considered smaller than xB
     */
    static bool isSmaller( const Seed& xA, const Seed& xB )
    {
        if( SeedLumping::getDelta( xA ) == SeedLumping::getDelta( xB ) )
            return xA.end( ) <= xB.start( );
        return SeedLumping::getDelta( xA ) < SeedLumping::getDelta( xB );
    }

    /**
     * @brief creates a new CompareSeedSets module.
     */
    CompareSeedSets( const ParameterSetManager& rParameters )
    {} // constructor

    std::shared_ptr<SeedSetComp> DLL_PORT( MA )
        execute( std::shared_ptr<Seeds> pGroundTruth, std::shared_ptr<Seeds> pData )
    {
        auto pRet = std::make_shared<SeedSetComp>( );
        auto xGroundTruthIt = pGroundTruth->begin( );
        auto xDataIt = pData->begin( );

        pRet->uiAmountGroundTruth += pGroundTruth->size( );
        pRet->uiAmountData += pData->size( );

        while( xGroundTruthIt != pGroundTruth->end( ) && xDataIt != pData->end( ) )
        {
            if( isSmaller( *xDataIt, *xGroundTruthIt ) )
                pRet->addData( xDataIt );
            else if( isSmaller( *xGroundTruthIt, *xDataIt ) )
                pRet->addGroundTruth( xGroundTruthIt );
            else
            {
                pRet->addOverlap( xDataIt, xGroundTruthIt );
                if( xDataIt->end( ) < xGroundTruthIt->end( ) )
                    pRet->addData( xDataIt );
                else
                    pRet->addGroundTruth( xGroundTruthIt );
            } // else
        } // while
        // remaining seeds in data
        while( xDataIt != pData->end( ) )
            pRet->addData( xDataIt );
        // remaining seeds in ground truth
        while( xGroundTruthIt != pGroundTruth->end( ) )
            pRet->addGroundTruth( xGroundTruthIt );

        return pRet;
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