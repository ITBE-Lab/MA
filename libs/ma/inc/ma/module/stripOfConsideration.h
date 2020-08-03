/**
 * @file stripOfConsideration.h
 * @brief Implements a Strip of Consideration.
 * @author Markus Schmidt
 */
#ifndef STRIP_OF_CONSIDERATION_H
#define STRIP_OF_CONSIDERATION_H

#include "ma/container/segment.h"
#include "ma/container/soc.h"
#include "ms/module/module.h"
#include <cmath>

namespace libMA
{

/**
 * @brief Used to quickly find areas with high density of @ref Seed "seeds".
 * @ingroup module
 */
class StripOfConsiderationSeeds : public libMS::Module<SoCPriorityQueue, false, Seeds, NucSeq, Pack, FMIndex>
{
  public:
    /**
     * @brief If the best SoC has seeds of accumulative length smaller than this, abort.
     * @details
     * Is multiplied by query length.
     * 0 = never abort.
     */
    const double fGiveUp;
    const size_t uiCurrHarmScoreMin;

    /**
     * @brief disable fGiveUp and fRelMinSeedSizeAmount if genome is too short
     */
    const nucSeqIndex uiMinGenomeSize;

    const size_t uiSoCWidth;

    inline static nucSeqIndex getPositionForBucketing( nucSeqIndex uiQueryLength, const Seed& xS )
    {
        return xS.start_ref( ) + ( uiQueryLength - xS.start( ) );
    } // function

    inline nucSeqIndex getStripSize( nucSeqIndex uiQueryLength, int iMatch, int iExtend, int iGap ) const
    {
        if( uiSoCWidth != 0 )
            return uiSoCWidth;

        return ( iMatch * uiQueryLength - iGap ) / iExtend;
    } // function


  public:
    StripOfConsiderationSeeds( const ParameterSetManager& rParameters )
        : fGiveUp( rParameters.getSelected( )->xHarmScoreMinRel->get( ) ),
          uiCurrHarmScoreMin( rParameters.getSelected( )->xHarmScoreMin->get( ) ),
          uiMinGenomeSize( rParameters.getSelected( )->xGenomeSizeDisable->get( ) ),
          uiSoCWidth( rParameters.getSelected( )->xSoCWidth->get( ) )
    {} // default constructor

    virtual std::shared_ptr<SoCPriorityQueue> DLL_PORT( MA ) execute( std::shared_ptr<Seeds> pSeeds,
                                                                      std::shared_ptr<NucSeq>
                                                                          pQuerySeq,
                                                                      std::shared_ptr<Pack>
                                                                          pRefSeq,
                                                                      std::shared_ptr<FMIndex>
                                                                          pFM_index );
}; // class


/**
 * @brief Filters a set of maximally extended seeds down to SMEMs
 * @ingroup module
 */
class ExtractSeeds : public libMS::Module<Seeds, false, SegmentVector, FMIndex, NucSeq, Pack>
{
    const unsigned int uiMaxAmbiguity;
    const size_t uiMinSeedSize;

    /**
     * @brief skip seeds with too much ambiguity
     * @details
     * True: skip all seeds with to much ambiguity
     * False: use max_hits instances of the seeds with more ambiguity
     */
    const bool bSkipLongBWTIntervals;

  public:
    /**
     * @brief extracts all seeds from a SegmentVector
     * @details
     * This mirrors seeds back to the reverse strand, so that the old harmonization & SoC code is not broken...
     */
    inline void emplaceAllNonBridgingSeed( SegmentVector& rSegmentVector, FMIndex& rxFM_index, Pack& rxRefSequence,
                                           Seeds& rvSeedVector, const nucSeqIndex uiQLen )
    {
        rSegmentVector.emplaceAllEachSeeds(
            rxFM_index, uiQLen, uiMaxAmbiguity, uiMinSeedSize, rvSeedVector,
            [&rxRefSequence, &rxFM_index, &rvSeedVector, &uiQLen]( ) {
                /*
                 * @note this bridging check is not required since we check weather a SoC
                 * is brigding in general.
                 * If any of the seeds within a SoC are bridging then the SoC is bridging.
                 * Could turn this into a debug assertion...
                 */
                rvSeedVector.back( ).uiDelta =
                    StripOfConsiderationSeeds::getPositionForBucketing( uiQLen, rvSeedVector.back( ) );

                size_t uiContig = rxRefSequence.uiSequenceIdForPosition( rvSeedVector.back( ).start_ref( ) );
                // move the delta position to the right by the contig index
                // this way we ensure seeds are ordered by contigs first
                // and then their delta positions within the contig.
                rvSeedVector.back( ).uiDelta += ( uiQLen + 10 ) * uiContig;
                // returning true since we want to continue extracting seeds
                return true;
            } // lambda
        );
    } // method
    ExtractSeeds( const ParameterSetManager& rParameters )
        : uiMaxAmbiguity( rParameters.getSelected( )->xMaximalSeedAmbiguity->get( ) ),
          uiMinSeedSize( rParameters.getSelected( )->xMinSeedLength->get( ) ),
          bSkipLongBWTIntervals( rParameters.getSelected( )->xSkipAmbiguousSeeds->get( ) )
    {} // constructor

    // overload
    virtual std::shared_ptr<Seeds> DLL_PORT( MA )
        execute( std::shared_ptr<SegmentVector> pSegments, std::shared_ptr<FMIndex> pFM_index,
                 std::shared_ptr<NucSeq> pQuerySeq, std::shared_ptr<Pack> pRefSeq )
    {

        const nucSeqIndex uiQLen = pQuerySeq->length( );

        // extract the seeds
        auto pSeeds = std::make_shared<Seeds>( );
        pSeeds->xStats.sName = pQuerySeq->sName;
        pSeeds->xStats.bSetMappingQualityToZero = pSegments->bSetMappingQualityToZero;
        // rough estimate of how many seeds we will have
        // (trying to avoid multiple allocations)
        pSeeds->reserve( pSegments->size( ) * 3 );

        emplaceAllNonBridgingSeed( *pSegments, // Segment vector (outcome of seeding)
                                   *pFM_index, *pRefSeq, *pSeeds, uiQLen ); // emplace all function call

        return pSeeds;
    } // method
}; // class

/**
 * @brief Used to quickly find areas with high density of @ref Seed "seeds".
 * @ingroup module
 */
class StripOfConsideration : public libMS::Module<SoCPriorityQueue, false, SegmentVector, NucSeq, Pack, FMIndex>
{
  public:
    StripOfConsiderationSeeds xHelper;
    ExtractSeeds xExtractHelper;

  public:
    StripOfConsideration( const ParameterSetManager& rParameters )
        : xHelper( rParameters ), xExtractHelper( rParameters )
    {} // default constructor

    virtual std::shared_ptr<SoCPriorityQueue> DLL_PORT( MA ) execute( std::shared_ptr<SegmentVector> pSegments,
                                                                      std::shared_ptr<NucSeq>
                                                                          pQuerySeq,
                                                                      std::shared_ptr<Pack>
                                                                          pRefSeq,
                                                                      std::shared_ptr<FMIndex>
                                                                          pFM_index );
}; // class


class GetAllFeasibleSoCs : public libMS::Module<Seeds, false, SoCPriorityQueue>
{
    const size_t uiSoCHeight;
    const nucSeqIndex uiMinNt;

  public:
    GetAllFeasibleSoCs( const ParameterSetManager& rParameters, nucSeqIndex uiMinNt )
        : uiSoCHeight( rParameters.getSelected( )->xSoCWidth->get( ) ), // same as width
          uiMinNt( uiMinNt )
    {} // constructor


    virtual std::shared_ptr<Seeds> DLL_PORT( MA ) execute( std::shared_ptr<SoCPriorityQueue> pSoCs )
    {
        auto pRet = std::make_shared<Seeds>( );

        while( !pSoCs->empty( ) && pSoCs->getScoreOfNextSoC( ) >= uiMinNt )
        {
            auto pSeeds = pSoCs->pop( );
            std::sort( pSeeds->begin( ), pSeeds->end( ),
                       []( auto& xA, auto& xB ) { return xA.start( ) < xB.start( ); } );
            nucSeqIndex uiNumNtLast = 0;
            nucSeqIndex uiMaxQ = 0;
            nucSeqIndex uiSizeLastSoC = 0;
            for( Seed& xSeed : *pSeeds )
            {
                if( xSeed.start( ) > uiMaxQ + uiSoCHeight )
                {
                    // last SoC had to little NT in it
                    if( uiNumNtLast < uiMinNt )
                        pRet->resize( pRet->size( ) - uiSizeLastSoC ); // remove it
                    uiNumNtLast = 0;
                    uiSizeLastSoC = 0;
                } // if
                uiNumNtLast += xSeed.size( );
                uiMaxQ = std::max( uiMaxQ, xSeed.end( ) );
                pRet->push_back( xSeed );
                uiSizeLastSoC++;
            } // for
            // very last SoC had to little NT in it
            if( uiNumNtLast < uiMinNt )
                pRet->resize( pRet->size( ) - uiSizeLastSoC ); // remove it
        } // while

        return pRet;
    } // method
}; // class


} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief export the bucketing @ref libMS::Module "module" to python.
 * @ingroup export
 */
#ifdef WITH_BOOST
void exportStripOfConsideration( );
#else
void exportStripOfConsideration( libMS::SubmoduleOrganizer& xOrganizer );
#endif
#endif

#endif