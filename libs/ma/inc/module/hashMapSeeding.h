/**
 * @file hashMapSeeding.h
 * @brief Computes the mapping quality of Alignments.
 * @author Markus Schmidt
 */
#pragma once

#include "container/alignment.h"
#include "module/harmonization.h"
#include "module/module.h"
#include <unordered_map>

namespace libMA
{
/**
 * @brief computes seeds between two sequences
 * @ingroup module
 * @details
 * Uses K-mer like seeding.
 * (this will produce many overlapping seeds)
 */
class HashMapSeeding : public libMS::Module<Seeds, false, NucSeq, NucSeq>
{
  public:
    nucSeqIndex uiSeedSize = 5;

    HashMapSeeding( const ParameterSetManager& rParameters )
    {} // constructor

    HashMapSeeding( )
    {} // constructor

    virtual nucSeqIndex minSeedSize( )
    {
        return uiSeedSize;
    } // method

    virtual nucSeqIndex getSeedSize( nucSeqIndex uiL1, nucSeqIndex uiL2 )
    {
        return minSeedSize( );
    } // method

    virtual std::unordered_multimap<std::string, size_t> getIndex( NucSeq& xQ1, nucSeqIndex uiSeedSize );
    virtual std::shared_ptr<Seeds> getSeeds( std::unordered_multimap<std::string, size_t> xIndex, NucSeq& xQ1,
                                             nucSeqIndex uiSeedSize );

    virtual std::shared_ptr<Seeds> EXPORTED execute( NucSeq& xQ1, NucSeq& xQ2 );
    virtual std::shared_ptr<Seeds> EXPORTED execute( std::shared_ptr<NucSeq> pQ1, std::shared_ptr<NucSeq> pQ2 )
    {
        return execute( *pQ1, *pQ2 );
    } // method

    virtual std::vector<std::shared_ptr<libMA::Seeds>> getAllSeeds( std::unordered_multimap<std::string, size_t> xIndex,
                                                                    std::vector<std::shared_ptr<NucSeq>> vIn )
    {
        std::vector<std::shared_ptr<libMA::Seeds>> vRet;
        for( auto pNucSeq : vIn )
            vRet.push_back( getSeeds( xIndex, *pNucSeq, xIndex.begin( )->first.size( ) ) );
        return vRet;
    } // method
}; // class

#if 0
/**
 * @brief fills in shorter seeds within gaps between seeds
 * @ingroup module
 * @details
 */
class ReSeeding : public libMS::Module<Seeds, false, Seeds, NucSeq, Pack>
{
  public:
    HashMapSeeding xHashMapSeeder;
    SeedLumping xLumper;
    nucSeqIndex uiPadding;

    ReSeeding( const ParameterSetManager& rParameters )
        : xHashMapSeeder( rParameters ),
          xLumper( rParameters ),
          uiPadding( rParameters.getSelected( )->xPaddingReSeeding->get( ) )
    {} // constructor

    virtual std::shared_ptr<Seeds> EXPORTED execute( std::shared_ptr<Seeds> pSeeds, std::shared_ptr<NucSeq> pQuery,
                                                     std::shared_ptr<Pack> pPack );

}; // class

/**
 * @brief fills seed sets using re seeding
 * @ingroup module
 * @details
 * Fills in shorter seeds within gaps.
 */
class FillSeedSet : public libMS::Module<Seeds, false, Seeds, NucSeq, FMIndex, Pack>
{
  public:
    HarmonizationSingle xSingle;
    ReSeeding xReseeding;

    FillSeedSet( const ParameterSetManager& rParameters ) : xSingle( rParameters ), xReseeding( rParameters )
    {} // default constructor

    // overload
    virtual std::shared_ptr<Seeds> EXPORTED execute( std::shared_ptr<Seeds> pSeedsIn,
                                                     std::shared_ptr<NucSeq>
                                                         pQuery,
                                                     std::shared_ptr<FMIndex>
                                                         pFMIndex,
                                                     std::shared_ptr<Pack>
                                                         pPack )
    {
        auto pAppend = std::make_shared<Seeds>( );
        // split seeds into forward and reverse strand
        auto pSecondaryStrand = pSeedsIn->extractStrand( false );
        // deal with seeds on forward strand
        while( !pSeedsIn->empty( ) )
            pAppend->append( xReseeding.execute( xSingle.execute( pSeedsIn, pQuery, pFMIndex ), pQuery, pPack ) );
        // deal with seeds on reverse strand
        while( !pSecondaryStrand->empty( ) )
            pAppend->append(
                xReseeding.execute( xSingle.execute( pSecondaryStrand, pQuery, pFMIndex ), pQuery, pPack ) );

        return pAppend;
    } // method
}; // class


/**
 * @brief extract seeds from a SoC priority queue
 * @ingroup module
 * @details
 * Extracts the seeds from a SoC priority queue,
 * then fills in shorter seeds within gaps.
 */
class ExtractFilledSeedSets
    : public libMS::Module<ContainerVector<std::shared_ptr<Seeds>>, false, SoCPriorityQueue, NucSeq, FMIndex, Pack>
{
  public:
    FillSeedSet xFill;

    /// @brief Extract at most x SoCs.
    const size_t uiMaxTries;

    ExtractFilledSeedSets( const ParameterSetManager& rParameters )
        : xFill( rParameters ), uiMaxTries( rParameters.getSelected( )->xMaxNumSoC->get( ) )
    {} // default constructor

    // overload
    virtual std::shared_ptr<ContainerVector<std::shared_ptr<Seeds>>>
        EXPORTED execute( std::shared_ptr<SoCPriorityQueue> pSoCsIn, std::shared_ptr<NucSeq> pQuery,
                          std::shared_ptr<FMIndex> pFMIndex, std::shared_ptr<Pack> pPack )
    {
        auto pSoCs = std::make_shared<ContainerVector<std::shared_ptr<Seeds>>>( );

        for( size_t uiNumTries = 0; uiNumTries < uiMaxTries && !pSoCsIn->empty( ); uiNumTries++ )
        {
            auto pAppend = xFill.execute( pSoCsIn->pop( ), pQuery, pFMIndex, pPack );

            if( !pAppend->empty( ) )
                pSoCs->push_back( pAppend );
        } // for

        return pSoCs;
    } // method
}; // class
#endif

} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief export the HashMapSeeding @ref libMA::Module "module" to python.
 * @ingroup export
 */
void exportHashMapSeeding( py::module& rxPyModuleId );
#endif
