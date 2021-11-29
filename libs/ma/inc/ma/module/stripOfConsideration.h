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
class StripOfConsiderationSeeds : public libMS::Module<SoCPriorityQueue, false, Seeds, NucSeq, Pack>
{
  public:
    /**
     * @brief If the best SoC has seeds of accumulative length smaller than this, abort.
     * @details
     * Is multiplied by query length.
     * 0 = never abort.
     */
    double fGiveUp;
    size_t uiCurrHarmScoreMin;

    /**
     * @brief disable fGiveUp and fRelMinSeedSizeAmount if genome is too short
     */
    const nucSeqIndex uiMinGenomeSize;

    const size_t uiSoCWidth;
    const bool bRectangular;

    inline static nucSeqIndex getPositionForBucketing( nucSeqIndex uiQueryLength, const Seed& xS, Pack& rxRefSequence,
                                                       bool bSplitStrands )
    {
        auto uiRet = xS.start_ref( );
        if( bSplitStrands && !xS.bOnForwStrand )
            // if strands are split reverse strand reads are sorted in reverse order so the position must be calculated
            // from the back (back of genome is 2*|R|+(|Q|+1)*num_contigs so that no seed's land in between contigs is
            // considered )
            uiRet = 2 * ( rxRefSequence.uiUnpackedSizeForwardStrand +
                          ( uiQueryLength + 1 ) * rxRefSequence.uiNumContigs( ) ) -
                    ( xS.start_ref( ) - xS.size( ) );
        return uiRet + ( uiQueryLength - xS.start( ) );
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
          uiSoCWidth( rParameters.getSelected( )->xSoCWidth->get( ) ),
          bRectangular( rParameters.getSelected( )->xRectangularSoc->get( ) )
    {} // default constructor

    virtual std::shared_ptr<SoCPriorityQueue> DLL_PORT( MA )
        execute( std::shared_ptr<Seeds> pSeeds, std::shared_ptr<NucSeq> pQuerySeq, std::shared_ptr<Pack> pRefSeq );
}; // class


/**
 * @brief Filters a set of maximally extended seeds down to SMEMs
 * @ingroup module
 */
class ExtractSeeds : public libMS::Module<Seeds, false, SegmentVector, FMIndex, NucSeq, Pack>
{
    const unsigned int uiMaxAmbiguity;
    const size_t uiMinSeedSize;
    const bool bRectangular;

    /**
     * @brief skip seeds with too much ambiguity
     * @details
     * True: skip all seeds with to much ambiguity
     * False: use max_hits instances of the seeds with more ambiguity
     */
    const bool bSkipLongBWTIntervals;

  public:
    static void setDeltaOfSeed( Seed& rSeed, const nucSeqIndex uiQLen, Pack& rxRefSequence, bool bSplitStrands )
    {
        rSeed.uiDelta =
            StripOfConsiderationSeeds::getPositionForBucketing( uiQLen, rSeed, rxRefSequence, bSplitStrands );

        size_t uiContig = rxRefSequence.uiSequenceIdForPosition( rSeed.start_ref( ) );

        if( bSplitStrands && !rSeed.bOnForwStrand )
            // if strands are split reverse strand reads are sorted in reverse order so the buffer space must be added
            // in reverse
            uiContig = rxRefSequence.uiNumContigs( ) - uiContig;
        // move the delta position to the right by the contig index
        // this way we ensure seeds are ordered by contigs first
        // and then their delta positions within the contig.
        rSeed.uiDelta += ( uiQLen + 1 ) * uiContig;
    } // method
    /**
     * @brief extracts all seeds from a SegmentVector
     * @details
     * This mirrors seeds back to the reverse strand, so that the old harmonization & SoC code is not broken...
     */
    inline void emplaceAllNonBridgingSeed( SegmentVector& rSegmentVector, FMIndex& rxFM_index, Pack& rxRefSequence,
                                           Seeds& rvSeedVector, const nucSeqIndex uiQLen )
    {
        rSegmentVector.emplaceAllEachSeeds( rxFM_index, uiQLen, uiMaxAmbiguity, uiMinSeedSize, rvSeedVector,
                                            [ & ]( ) {
                                                setDeltaOfSeed( rvSeedVector.back( ), uiQLen, rxRefSequence,
                                                                !bRectangular );
                                                // returning true since we want to continue extracting seeds
                                                return true;
                                            } // lambda
        );
    } // method
    ExtractSeeds( const ParameterSetManager& rParameters )
        : uiMaxAmbiguity( rParameters.getSelected( )->xMaximalSeedAmbiguity->get( ) ),
          uiMinSeedSize( rParameters.getSelected( )->xMinSeedLength->get( ) ),
          bRectangular( rParameters.getSelected( )->xRectangularSoc->get( ) ),
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

  public:
    const nucSeqIndex uiMinNt;
    const double dMinNtRelative;

    GetAllFeasibleSoCs( const ParameterSetManager& rParameters )
        : uiSoCHeight( rParameters.getSelected( )->xSoCWidth->get( ) ), // same as width
          uiMinNt( rParameters.getSelected( )->xMinNtInSoc->get( ) ),
          dMinNtRelative( rParameters.getSelected( )->xMinNtInSocRelative->get( ) )
    {} // constructor


    virtual std::shared_ptr<Seeds> DLL_PORT( MA ) execute( std::shared_ptr<SoCPriorityQueue> pSoCs )
    {
        assert( pSoCs->uiQLen > 0 );
        const nucSeqIndex uiThreshold = std::min(uiMinNt, (nucSeqIndex) (dMinNtRelative*pSoCs->uiQLen));
        auto pRet = std::make_shared<Seeds>( );

        while( !pSoCs->empty( ) && pSoCs->getScoreOfNextSoC( ) >= uiThreshold )
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
                    if( uiNumNtLast < uiThreshold )
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
            if( uiNumNtLast < uiThreshold )
                pRet->resize( pRet->size( ) - uiSizeLastSoC ); // remove it
        } // while

        return pRet;
    } // method
}; // class

class GetAllFeasibleSoCsAsSet : public libMS::Module<SeedsSet, false, SoCPriorityQueue>
{
    const size_t uiSoCHeight;

  public:
    nucSeqIndex uiMinNt;
    const double dMinNtRelative;

    GetAllFeasibleSoCsAsSet( const ParameterSetManager& rParameters )
        : uiSoCHeight( rParameters.getSelected( )->xSoCWidth->get( ) ), // same as width
          uiMinNt( rParameters.getSelected( )->xMinNtInSoc->get( ) ),
          dMinNtRelative( rParameters.getSelected( )->xMinNtInSocRelative->get( ) )
    {} // constructor


    virtual std::shared_ptr<SeedsSet> DLL_PORT( MA ) execute( std::shared_ptr<SoCPriorityQueue> pSoCs )
    {
        assert( pSoCs->uiQLen > 0 );
        const nucSeqIndex uiThreshold = std::min(uiMinNt, (nucSeqIndex) (dMinNtRelative*pSoCs->uiQLen));
        auto pRet = std::make_shared<SeedsSet>( );

        while( !pSoCs->empty( ) && pSoCs->getScoreOfNextSoC( ) >= uiThreshold )
        {
#if 0 // turn on/off the splitting of SoCs on gaps
            auto pSeeds = pSoCs->pop( );
            pRet->xContent.push_back( pSeeds );
#else
            pRet->xContent.push_back( std::make_shared<Seeds>( ) );
            auto pSeeds = pSoCs->pop( );
            std::sort( pSeeds->begin( ), pSeeds->end( ),
                       []( auto& xA, auto& xB ) { return xA.start( ) < xB.start( ); } );
            nucSeqIndex uiNumNtLast = 0;
            nucSeqIndex uiMaxQ = pSeeds->front( ).end( );
            for( Seed& xSeed : *pSeeds )
            {
                if( xSeed.start( ) > uiMaxQ + uiSoCHeight )
                {
                    // last SoC had to little NT in it
                    if( uiNumNtLast < uiThreshold )
                        pRet->xContent.pop_back( ); // remove it
                    pRet->xContent.push_back( std::make_shared<Seeds>( ) );
                    uiNumNtLast = 0;
                } // if
                uiNumNtLast += xSeed.size( );
                uiMaxQ = std::max( uiMaxQ, xSeed.end( ) );
                pRet->xContent.back( )->push_back( xSeed );
            } // for
            // very last SoC had to little NT in it
            if( uiNumNtLast < uiThreshold )
                pRet->xContent.pop_back( ); // remove it
#endif
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