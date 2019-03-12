/**
 * @file hashMapSeeding.h
 * @brief Computes the mapping quality of Alignments.
 * @author Markus Schmidt
 */
#pragma once

#include "container/alignment.h"
#include "module/harmonization.h"
#include "module/module.h"

namespace libMA
{
/**
 * @brief computes seeds between two sequences
 * @ingroup module
 * @details
 * Uses K-mer like seeding.
 * (this will produce many overlapping seeds)
 */
class HashMapSeeding : public Module<Seeds, false, NucSeq, NucSeq>
{
  public:
    const size_t uiSeedSize;

    HashMapSeeding( const ParameterSetManager& rParameters )
        : uiSeedSize( rParameters.getSelected( )->xSecSeedSize->get( ) )
    {} // constructor

    virtual std::shared_ptr<Seeds> EXPORTED execute( std::shared_ptr<NucSeq> pQ1, std::shared_ptr<NucSeq> pQ2 );

}; // class

/**
 * @brief fills in shorter seeds within gaps between seeds
 * @ingroup module
 * @details
 */
class ReSeeding : public Module<Seeds, false, Seeds, NucSeq, Pack>
{
  public:
    HashMapSeeding xHashMapSeeder;
    SeedLumping xLumper;

    ReSeeding( const ParameterSetManager& rParameters ) : xHashMapSeeder( rParameters ), xLumper( rParameters )
    {} // constructor

    virtual std::shared_ptr<Seeds> EXPORTED execute( std::shared_ptr<Seeds> pSeeds, std::shared_ptr<NucSeq> pQuery,
                                                     std::shared_ptr<Pack> pPack );

}; // class


/**
 * @brief extract seeds from a SoC priority queue
 * @ingroup module
 * @details
 * Extracts the seeds from a SoC priority queue,
 * then fills in shorter seeds within gaps.
 */
class ExtractFilledSeedSets
    : public Module<ContainerVector<std::shared_ptr<Seeds>>, false, SoCPriorityQueue, NucSeq, FMIndex, Pack>
{
  public:
    HarmonizationSingle xSingle;
    ReSeeding xReseeding;

    /// @brief Extract at most x SoCs.
    const size_t uiMaxTries;

    ExtractFilledSeedSets( const ParameterSetManager& rParameters )
        : xSingle( rParameters ),
          xReseeding( rParameters ),
          uiMaxTries( rParameters.getSelected( )->xMaxNumSoC->get( ) )
    {} // default constructor

    // overload
    virtual std::shared_ptr<ContainerVector<std::shared_ptr<Seeds>>>
        EXPORTED execute( std::shared_ptr<SoCPriorityQueue> pSoCsIn, std::shared_ptr<NucSeq> pQuery,
                          std::shared_ptr<FMIndex> pFMIndex, std::shared_ptr<Pack> pPack )
    {
        auto pSoCs = std::make_shared<ContainerVector<std::shared_ptr<Seeds>>>( );

        for( size_t uiNumTries = 0; uiNumTries < uiMaxTries && !pSoCsIn->empty( ); uiNumTries++ )
        {
            auto pAppend = std::make_shared<Seeds>( );
            auto pSoC = pSoCsIn->pop( );
            // split seeds into forward and reverse strand
            auto pSecondaryStrand = pSoC->extractStrand( false );
            // deal with seeds on forward strand
            while( !pSoC->empty( ) )
                pAppend->append( xReseeding.execute( xSingle.execute( pSoC, pQuery, pFMIndex ), pQuery, pPack ) );
            // deal with seeds on reverse strand
            while( !pSecondaryStrand->empty( ) )
                pAppend->append(
                    xReseeding.execute( xSingle.execute( pSecondaryStrand, pQuery, pFMIndex ), pQuery, pPack ) );

            if( !pAppend->empty( ) )
                pSoCs->push_back( pAppend );
        } // for

        return pSoCs;
    } // method
}; // class

} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief export the HashMapSeeding @ref Module "module" to python.
 * @ingroup export
 */
void exportHashMapSeeding( py::module& rxPyModuleId );
#endif
